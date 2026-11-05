"""Static oil-and-gas basin boxes used to choose Google News RSS queries."""
from typing import List, NamedTuple, Sequence

from backend.schemas import BBox


class Basin(NamedTuple):
    id: str
    name: str
    bbox: BBox
    aliases: Sequence[str]


BASINS: Sequence[Basin] = (
    Basin("permian", "Permian", BBox(-105.5, 30.0, -100.0, 33.8), ("Permian", "Delaware Basin")),
    Basin("eagle_ford", "Eagle Ford", BBox(-101.5, 27.4, -96.5, 29.9), ("Eagle Ford",)),
    Basin("haynesville", "Haynesville", BBox(-96.0, 31.0, -92.5, 33.2), ("Haynesville",)),
    Basin("marcellus", "Marcellus", BBox(-82.5, 37.5, -74.5, 42.8), ("Marcellus",)),
    Basin("utica", "Utica", BBox(-83.0, 39.0, -78.0, 42.6), ("Utica",)),
    Basin("bakken", "Bakken", BBox(-105.5, 46.5, -100.0, 49.2), ("Bakken",)),
    Basin("dj", "DJ", BBox(-105.5, 39.5, -102.0, 41.6), ("DJ", "Denver-Julesburg", "Denver Julesburg")),
    Basin("anadarko", "Anadarko", BBox(-102.5, 34.5, -96.5, 37.6), ("Anadarko",)),
    Basin("uinta", "Uinta", BBox(-111.5, 39.0, -108.5, 41.0), ("Uinta",)),
    Basin("montney", "Montney", BBox(-123.5, 54.0, -118.0, 57.5), ("Montney",)),
    Basin("duvernay", "Duvernay", BBox(-118.0, 52.0, -112.0, 55.5), ("Duvernay",)),
    Basin("wcsb", "WCSB", BBox(-125.0, 48.5, -108.0, 60.0), ("WCSB", "Western Canada Sedimentary")),
    Basin("appalachian", "Appalachian", BBox(-83.0, 37.0, -74.0, 42.8), ("Appalachian",)),
)


def bbox_overlaps(left: BBox, right: BBox) -> bool:
    return left.west <= right.east and left.east >= right.west and left.south <= right.north and left.north >= right.south


def overlapping(bbox: BBox) -> List[Basin]:
    return [basin for basin in BASINS if bbox_overlaps(bbox, basin.bbox)]
