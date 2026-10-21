import pytest
from ttm.spikeinterface import match_stream

IMU_RIG = (
    [
        "Record Node 102#Acquisition_Board-101.IMU_port_A",
        "Record Node 102#Acquisition_Board-101.acquisition_board",
        "Record Node 102#Acquisition_Board-101.memory_usage",
    ],
    ["0", "1", "2"],
)
TTM1_RIG = (
    [
        "Record Node 121#Acquisition_Board-120.acquisition_board",
        "Record Node 121#Acquisition_Board-120.memory_usage",
    ],
    ["0", "1"],
)


def test_suffix_match_is_rig_independent() -> None:
    assert match_stream(*IMU_RIG, "acquisition_board") == "1"
    assert match_stream(*TTM1_RIG, "acquisition_board") == "0"
    assert match_stream(*IMU_RIG, "IMU_port_A") == "0"


def test_full_name_and_no_partial_suffix() -> None:
    assert (
        match_stream(*IMU_RIG, "Record Node 102#Acquisition_Board-101.memory_usage")
        == "2"
    )
    with pytest.raises(ValueError, match="matches 0 streams"):
        match_stream(*IMU_RIG, "board")


def test_ambiguous_match_raises() -> None:
    names = [
        "Record Node 1#A-100.acquisition_board",
        "Record Node 1#B-200.acquisition_board",
    ]
    with pytest.raises(ValueError, match="matches 2 streams"):
        match_stream(names, ["0", "1"], "acquisition_board")
    assert match_stream(names, ["0", "1"], names[1]) == "1"
