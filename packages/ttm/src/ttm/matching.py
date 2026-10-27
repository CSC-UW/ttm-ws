"""Template matching run on each channel group (tetrode) on its own.

SpikeInterface's matchers (circus-omp, wobble) fit every template to every
channel of a chunk together: the candidate spikes of a pass, and the decision
to stop, are shared by the whole probe. ``PerGroupMatching`` reads a chunk once
and runs one matcher per group on that group's channels with that group's
templates, so each tetrode gets its own candidates and its own stopping. It
needs group-sparse templates (each template's channels inside one group).
``find_spikes_per_group`` has the signature of
``spikeinterface.sortingcomponents.matching.find_spikes_from_templates``.
"""

from __future__ import annotations

import numpy as np
from spikeinterface.core import Templates
from spikeinterface.core.node_pipeline import run_node_pipeline
from spikeinterface.sortingcomponents.matching.base import (
    BaseTemplateMatching,
    _base_matching_dtype,
)
from spikeinterface.sortingcomponents.matching.method_list import matching_methods


def template_groups(templates: Templates, channel_groups: np.ndarray) -> np.ndarray:
    """Group of each template, from its sparsity mask; every template must sit inside one group."""
    if templates.sparsity_mask is None:
        raise ValueError("per-group matching needs group-sparse templates")
    out = np.empty(templates.num_units, dtype=channel_groups.dtype)
    for i, m in enumerate(templates.sparsity_mask):
        g = np.unique(channel_groups[m])
        if len(g) != 1:
            raise ValueError(f"template {templates.unit_ids[i]} spans groups {g.tolist()}")
        out[i] = g[0]
    return out


def group_templates(templates: Templates, units: np.ndarray, channels: np.ndarray) -> Templates:
    """The templates ``units`` restricted to ``channels`` (indices), sparse with the restricted mask."""
    dense = templates.get_dense_templates()[units][:, :, channels].astype(np.float32)
    mask = templates.sparsity_mask[units][:, channels]
    packed = np.zeros((len(units), dense.shape[1], int(mask.sum(axis=1).max())), dtype=np.float32)
    for i, m in enumerate(mask):
        packed[i, :, : m.sum()] = dense[i][:, m]
    return Templates(templates_array=packed, sampling_frequency=templates.sampling_frequency,
                     nbefore=templates.nbefore, is_in_uV=templates.is_in_uV, sparsity_mask=mask,
                     channel_ids=np.asarray(templates.channel_ids)[channels],
                     unit_ids=np.asarray(templates.unit_ids)[units])


class PerGroupMatching(BaseTemplateMatching):
    """One ``method`` matcher per channel group, run on the same chunk; spikes in the full bank's indices."""

    def __init__(self, recording, templates, return_output=True, method="circus-omp", group_property="group",
                 **method_kwargs):
        BaseTemplateMatching.__init__(self, recording, templates, return_output=return_output)
        groups = np.asarray(recording.get_property(group_property))
        if not np.array_equal(np.asarray(templates.channel_ids), np.asarray(recording.channel_ids)):
            raise ValueError("templates and recording must have the same channels in the same order")
        unit_group = template_groups(templates, groups)
        self.parts = []
        for g in np.unique(unit_group):
            channels = np.flatnonzero(groups == g)
            units = np.flatnonzero(unit_group == g)
            sub = recording.select_channels(np.asarray(recording.channel_ids)[channels])
            node = matching_methods[method](sub, templates=group_templates(templates, units, channels),
                                            **method_kwargs)
            self.parts.append((channels, units, node))
        self.margin = max(node.get_margin() for _, _, node in self.parts)

    def get_margin(self):
        return self.margin

    def compute_matching(self, traces, start_frame, end_frame, segment_index):
        found = []
        for channels, units, node in self.parts:
            s = node.compute_matching(np.ascontiguousarray(traces[:, channels]), start_frame, end_frame,
                                      segment_index)
            s["cluster_index"] = units[s["cluster_index"]]
            s["channel_index"] = channels[s["channel_index"]]
            found.append(s)
        spikes = np.concatenate(found)
        return spikes[np.argsort(spikes["sample_index"], kind="stable")]

    def clean(self):
        for _, _, node in self.parts:
            node.clean()


def find_spikes_per_group(recording, templates, method=None, method_kwargs=None, extra_outputs=False,
                          pipeline_kwargs=None, verbose=False, job_kwargs=None, group_property="group"):
    """``find_spikes_from_templates`` with each channel group matched on its own (``PerGroupMatching``)."""
    method_kwargs = dict(method_kwargs or {})
    method = method_kwargs.pop("method", method)
    if extra_outputs:
        raise NotImplementedError("extra_outputs is not available per group")
    if len(templates.unit_ids) == 0:
        return np.zeros(0, dtype=_base_matching_dtype)
    node = PerGroupMatching(recording, templates, method=method, group_property=group_property, **method_kwargs)
    spikes = run_node_pipeline(recording, [node], job_kwargs or {}, job_name=f"find spikes per group ({method})",
                               squeeze_output=True, names=["spikes"], verbose=verbose, **(pipeline_kwargs or {}))
    node.clean()
    return spikes
