import tempfile
import unittest
from datetime import date

from backend.schemas import BBox, DateRange
from backend.sources import ghgrp_ca, ghgrp_us
from backend.sources.file_cache import co2e_tonnes_to_ch4_kg_hr, tonnes_per_year_to_kg_hr

US_CSV = (
    "facility_id,latitude,longitude,city,state,facility_name,year,reported_subparts,"
    "parent_company,sub_part_id,co2e_emission,gas_code\n"
    "1000001,48.828707,-122.685533,FERNDALE,WA,PSE Ferndale Generating Station,2023,D,"
    "PUGET HOLDINGS LLC,2,219.15,CH4\n"
    "1000001,48.828707,-122.685533,FERNDALE,WA,PSE Ferndale Generating Station,2023,D,"
    "PUGET HOLDINGS LLC,1,219.15,CH4\n"
    "2000002,51.04,-114.07,CALGARY,AB,Outside,2023,HH,TEST,1,1000,CH4\n"
)
US_BBOX = BBox(west=-123.0, south=48.0, east=-122.0, north=49.0)

CA_CSV = (
    "GHGRP ID No. / No d'identification du PDGES,Reference Year / Année de référence,"
    "Facility Name / Nom de l'installation,Latitude,Longitude,"
    "English Facility NAICS Code Description / Description du code SCIAN de l'installation en anglais,"
    "CH4 (tonnes),CH4 (tonnes CO2e / tonnes éq. CO2)\n"
    "G1,2024,Keele Valley Landfill,43.85,-79.50,Waste treatment and disposal,87.66,2191.5\n"
    "G1,2023,Keele Valley Landfill,43.85,-79.50,Waste treatment and disposal,10.0,250.0\n"
    "G2,2024,Outside Mine,51.04,-114.07,Oil and gas extraction,5.0,125.0\n"
    "G3,2024,Zero CH4,43.90,-79.40,Electric power,0,0\n"
)
CA_BBOX = BBox(west=-80.0, south=43.0, east=-79.0, north=44.0)
DATES = DateRange(date(2026, 1, 1), date(2026, 1, 31))


class GhgrpUsTests(unittest.TestCase):
    def test_aggregates_subparts_and_converts_co2e(self):
        plumes = ghgrp_us.parse_csv(US_CSV)
        self.assertEqual(len(plumes), 2)
        ferndale = next(plume for plume in plumes if plume.id == "1000001")
        self.assertEqual(ferndale.source, "ghgrp_us")
        self.assertEqual(ferndale.observed_at.year, 2023)
        self.assertEqual(ferndale.instrument, "US GHGRP")
        self.assertAlmostEqual(ferndale.emission_kg_hr, co2e_tonnes_to_ch4_kg_hr(438.3))
        self.assertEqual(ferndale.details["GWP"], "25")
        self.assertEqual(ferndale.details["Name"], "PSE Ferndale Generating Station")

    def test_fetch_filters_bbox(self):
        with tempfile.TemporaryDirectory() as tmp:
            collection = ghgrp_us.fetch_plumes(
                US_BBOX,
                DATES,
                cache_dir=tmp,
                fetch=lambda url: US_CSV.encode("utf-8"),
            )
        self.assertEqual([plume.id for plume in collection.items], ["1000001"])


class GhgrpCaTests(unittest.TestCase):
    def test_keeps_latest_year_with_ch4(self):
        plumes = ghgrp_ca.parse_csv(CA_CSV)
        ids = [plume.id for plume in plumes]
        self.assertEqual(ids, ["G1", "G2"])
        landfill = plumes[0]
        self.assertEqual(landfill.source, "ghgrp_ca")
        self.assertEqual(landfill.observed_at.year, 2024)
        self.assertAlmostEqual(landfill.emission_kg_hr, tonnes_per_year_to_kg_hr(87.66))
        self.assertEqual(landfill.sector, "Waste treatment and disposal")
        self.assertEqual(landfill.instrument, "Canada GHGRP")

    def test_fetch_filters_bbox(self):
        with tempfile.TemporaryDirectory() as tmp:
            collection = ghgrp_ca.fetch_plumes(
                CA_BBOX,
                DATES,
                cache_dir=tmp,
                fetch=lambda url: CA_CSV.encode("utf-8"),
            )
        self.assertEqual([plume.id for plume in collection.items], ["G1"])


if __name__ == "__main__":
    unittest.main()
