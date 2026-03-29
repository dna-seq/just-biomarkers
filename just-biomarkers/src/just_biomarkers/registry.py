"""Clock registry -- static definitions of all supported methylation clocks.

Each entry mirrors the biolearn ``model_definitions`` dict but is stored as
a :class:`ClockDefinition` pydantic model so it can be serialised, filtered,
and introspected without importing heavy dependencies.
"""
from __future__ import annotations

from just_biomarkers.models import ClockDefinition

CLOCK_DEFINITIONS: dict[str, ClockDefinition] = {
    "Horvathv1": ClockDefinition(
        name="Horvathv1",
        year=2013,
        tissue="Multi-tissue",
        source="https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115",
        output="Age (Years)",
        coefficient_file="Horvath1.csv",
        transform_name="horvathv1",
    ),
    "Hannum": ClockDefinition(
        name="Hannum",
        year=2013,
        tissue="Blood",
        source="https://www.sciencedirect.com/science/article/pii/S1097276512008933",
        output="Age (Years)",
        coefficient_file="Hannum.csv",
    ),
    "PhenoAge": ClockDefinition(
        name="PhenoAge",
        year=2018,
        tissue="Blood",
        source="https://www.aging-us.com/article/101414/text",
        output="Age (Years)",
        coefficient_file="PhenoAge.csv",
    ),
    "Lin": ClockDefinition(
        name="Lin",
        year=2016,
        tissue="Blood",
        source="https://www.aging-us.com/article/100908/text",
        output="Age (Years)",
        coefficient_file="Lin.csv",
    ),
    "Horvathv2": ClockDefinition(
        name="Horvathv2",
        year=2018,
        tissue="Skin + blood",
        source="https://www.aging-us.com/article/101508/text",
        output="Age (Years)",
        coefficient_file="Horvath2.csv",
        transform_name="horvathv2",
    ),
    "DunedinPoAm38": ClockDefinition(
        name="DunedinPoAm38",
        year=2020,
        tissue="Blood",
        source="https://elifesciences.org/articles/54870#s2",
        output="Aging Rate (Years/Year)",
        coefficient_file="DunedinPoAm38.csv",
    ),
    "DunedinPACE": ClockDefinition(
        name="DunedinPACE",
        year=2022,
        tissue="Blood",
        source="https://www.proquest.com/docview/2634411178",
        output="Aging Rate (Years/Year)",
        coefficient_file="DunedinPACE.csv",
        preprocess_name="dunedin_pace",
        default_imputation="none",
    ),
    "Zhang_10": ClockDefinition(
        name="Zhang_10",
        year=2019,
        tissue="Blood",
        source="https://www.nature.com/articles/ncomms14617",
        output="Mortality Risk",
        coefficient_file="Zhang_10.csv",
    ),
    "YingCausAge": ClockDefinition(
        name="YingCausAge",
        year=2022,
        tissue="Blood",
        source="https://www.biorxiv.org/content/10.1101/2022.10.07.511382v2",
        output="Age (Years)",
        coefficient_file="YingCausAge.csv",
    ),
    "YingDamAge": ClockDefinition(
        name="YingDamAge",
        year=2022,
        tissue="Blood",
        source="https://www.biorxiv.org/content/10.1101/2022.10.07.511382v2",
        output="Age (Years)",
        coefficient_file="YingDamAge.csv",
    ),
    "YingAdaptAge": ClockDefinition(
        name="YingAdaptAge",
        year=2022,
        tissue="Blood",
        source="https://www.biorxiv.org/content/10.1101/2022.10.07.511382v2",
        output="Age (Years)",
        coefficient_file="YingAdaptAge.csv",
    ),
    "PEDBE": ClockDefinition(
        name="PEDBE",
        year=2019,
        tissue="Buccal",
        source="https://www.pnas.org/doi/10.1073/pnas.1820843116",
        output="Age (Years)",
        coefficient_file="PEDBE.csv",
        transform_name="pedbe",
    ),
    "DNAmClockCortical": ClockDefinition(
        name="DNAmClockCortical",
        year=2020,
        tissue="Human Cortex",
        source="https://doi.org/10.1093/brain/awaa334",
        output="Human Cortex Age (Years)",
        coefficient_file="DNAmClockCortical.csv",
        transform_name="cortical",
    ),
    "VidalBralo": ClockDefinition(
        name="VidalBralo",
        year=2018,
        tissue="Blood",
        source="https://doi.org/10.3389/fgene.2016.00126",
        output="Age (Years)",
        coefficient_file="VidalBralo.csv",
        transform_name="vidal_bralo",
    ),
    "DNAmTL": ClockDefinition(
        name="DNAmTL",
        year=2019,
        tissue="Blood, Adipose",
        source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6738410/",
        output="Telomere Length",
        coefficient_file="DNAmTL.csv",
    ),
    "AlcoholMcCartney": ClockDefinition(
        name="AlcoholMcCartney",
        year=2018,
        tissue="Blood",
        source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6158884/",
        output="Alcohol Consumption",
        coefficient_file="Alcohol.csv",
    ),
    "SmokingMcCartney": ClockDefinition(
        name="SmokingMcCartney",
        year=2018,
        tissue="Blood",
        source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6158884/",
        output="Smoking Status",
        coefficient_file="Smoking.csv",
    ),
    "BMI_McCartney": ClockDefinition(
        name="BMI_McCartney",
        year=2018,
        tissue="Blood",
        source="https://doi.org/10.1186/s13059-018-1514-1",
        output="BMI",
        coefficient_file="BMI_McCartney.csv",
        transform_name="sigmoid",
    ),
    "EducationMcCartney": ClockDefinition(
        name="EducationMcCartney",
        year=2018,
        tissue="Blood",
        source="https://doi.org/10.1186/s13059-018-1514-1",
        output="Educational Attainment",
        coefficient_file="EducationMcCartney.csv",
        transform_name="sigmoid",
    ),
    "BodyFatMcCartney": ClockDefinition(
        name="BodyFatMcCartney",
        year=2018,
        tissue="Blood",
        source="https://doi.org/10.1186/s13059-018-1514-1",
        output="Percentage Body Fat",
        coefficient_file="BodyFatMcCartney.csv",
        transform_name="sigmoid",
    ),
    "HDLCholesterolMcCartney": ClockDefinition(
        name="HDLCholesterolMcCartney",
        year=2018,
        tissue="Blood",
        source="https://doi.org/10.1186/s13059-018-1514-1",
        output="HDL Cholesterol",
        coefficient_file="HDLCholesterolMcCartney.csv",
        transform_name="sigmoid",
    ),
    "LDLCholesterolMcCartney": ClockDefinition(
        name="LDLCholesterolMcCartney",
        year=2018,
        tissue="Blood",
        source="https://doi.org/10.1186/s13059-018-1514-1",
        output="LDL with Remnant Cholesterol",
        coefficient_file="LDLCholesterolMcCartney.csv",
        transform_name="sigmoid",
    ),
    "TotalCholesterolMcCartney": ClockDefinition(
        name="TotalCholesterolMcCartney",
        year=2018,
        tissue="Blood",
        source="https://doi.org/10.1186/s13059-018-1514-1",
        output="Total Cholesterol",
        coefficient_file="TotalCholesterolMcCartney.csv",
        transform_name="sigmoid",
    ),
    "StocZ": ClockDefinition(
        name="StocZ",
        year=2024,
        tissue="Blood",
        source="https://doi.org/10.1038/s43587-024-00600-8",
        output="Mortality Risk",
        coefficient_file="StocZ.csv",
        transform_name="stocz",
    ),
    "StocP": ClockDefinition(
        name="StocP",
        year=2024,
        tissue="Blood",
        source="https://doi.org/10.1038/s43587-024-00600-8",
        output="Age (Years)",
        coefficient_file="StocP.csv",
        transform_name="stocp",
    ),
    "StocH": ClockDefinition(
        name="StocH",
        year=2024,
        tissue="Multi-tissue",
        source="https://doi.org/10.1038/s43587-024-00600-8",
        output="Age (Years)",
        coefficient_file="StocH.csv",
        transform_name="stoch",
    ),
    "HRSInCHPhenoAge": ClockDefinition(
        name="HRSInCHPhenoAge",
        year=2022,
        tissue="Blood",
        source="https://www.nature.com/articles/s43587-022-00248-2",
        output="Age (Years)",
        coefficient_file="HRSInCHPhenoAge.csv",
    ),
    "EpiTOC1": ClockDefinition(
        name="EpiTOC1",
        year=2016,
        tissue="Blood",
        source="https://doi.org/10.1186/s13059-016-1064-3",
        output="Stem Cell Division Rate",
        coefficient_file="EpiTOC1.csv",
    ),
    "Knight": ClockDefinition(
        name="Knight",
        year=2016,
        tissue="Cord Blood",
        source="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1068-z",
        output="Gestational Age",
        coefficient_file="Knight.csv",
    ),
    "LeeControl": ClockDefinition(
        name="LeeControl",
        year=2019,
        tissue="Placenta",
        source="https://www.aging-us.com/article/102049/text",
        output="Gestational Age",
        coefficient_file="LeeControl.csv",
    ),
    "LeeRefinedRobust": ClockDefinition(
        name="LeeRefinedRobust",
        year=2019,
        tissue="Placenta",
        source="https://www.aging-us.com/article/102049/text",
        output="Gestational Age",
        coefficient_file="LeeRefinedRobust.csv",
    ),
    "LeeRobust": ClockDefinition(
        name="LeeRobust",
        year=2019,
        tissue="Placenta",
        source="https://www.aging-us.com/article/102049/text",
        output="Gestational Age",
        coefficient_file="LeeRobust.csv",
    ),
    "BMI_Reed": ClockDefinition(
        name="BMI_Reed",
        year=2020,
        tissue="Blood",
        source="https://doi.org/10.1186/s13148-020-00841-5",
        output="BMI",
        coefficient_file="BMI_Reed.csv",
    ),
    "CVD_Westerman": ClockDefinition(
        name="CVD_Westerman",
        year=2020,
        tissue="Blood",
        source="https://doi.org/10.1161/JAHA.119.015299",
        output="Coronary Heart Disease Status",
        coefficient_file="CVD_Westermann.csv",
        transform_name="sigmoid",
    ),
    "AD_Bahado-Singh": ClockDefinition(
        name="AD_Bahado-Singh",
        year=2021,
        tissue="Blood",
        source="https://doi.org/10.1371/journal.pone.0248375",
        output="Alzheimer's Disease Status",
        coefficient_file="AD_Bahado-Singh.csv",
        transform_name="ad_bahado_singh",
    ),
    "DepressionBarbu": ClockDefinition(
        name="DepressionBarbu",
        year=2021,
        tissue="Blood",
        source="https://doi.org/10.1038/s41380-020-0808-3",
        output="Depression Risk",
        coefficient_file="DepressionBarbu.csv",
    ),
    "Weidner": ClockDefinition(
        name="Weidner",
        year=2014,
        tissue="Blood",
        source="https://doi.org/10.1186/gb-2014-15-2-r24",
        output="Age (Years)",
        coefficient_file="Weidner.csv",
        transform_name="weidner",
    ),
    "Garagnani": ClockDefinition(
        name="Garagnani",
        year=2012,
        tissue="Blood",
        source="https://pubmed.ncbi.nlm.nih.gov/23061750/",
        output="Age (Years)",
        coefficient_file="Garagnani.csv",
    ),
    "Mayne": ClockDefinition(
        name="Mayne",
        year=2016,
        tissue="Placenta",
        source="https://doi.org/10.2217/epi-2016-0103",
        output="Gestational Age",
        coefficient_file="Mayne.csv",
        transform_name="mayne",
    ),
    "ProstateCancerKirby": ClockDefinition(
        name="ProstateCancerKirby",
        year=2017,
        tissue="Prostate",
        source="https://doi.org/10.1186/s12885-017-3252-2",
        output="Prostate Cancer Status",
        coefficient_file="ProstateCancerKirby.csv",
    ),
    "HepatoXu": ClockDefinition(
        name="HepatoXu",
        year=2017,
        tissue="Circulating DNA",
        source="https://doi.org/10.1038/nmat4997",
        output="Hepatocellular Carcinoma Status",
        coefficient_file="HepatoXu.csv",
    ),
    "Bohlin": ClockDefinition(
        name="Bohlin",
        year=2017,
        tissue="Cord Blood",
        source="https://doi.org/10.1186/s13059-016-1063-4",
        output="Age (days)",
        coefficient_file="Bohlin.csv",
        transform_name="bohlin",
    ),
    "Bocklandt": ClockDefinition(
        name="Bocklandt",
        year=2011,
        tissue="Blood",
        source="https://doi.org/10.1371/journal.pone.0014821",
        output="Age (Years)",
        coefficient_file="Bocklandt.csv",
    ),
    "DownSyndrome": ClockDefinition(
        name="DownSyndrome",
        year=2021,
        tissue="Blood",
        source="https://www.nature.com/articles/s41467-021-21064-z",
        output="Down Syndrome Prediction",
        coefficient_file="down_syndrome.csv",
        default_imputation="averaging",
    ),
}


def list_clocks() -> list[str]:
    """Return sorted names of all registered clocks."""
    return sorted(CLOCK_DEFINITIONS.keys())


def get_clock(name: str) -> ClockDefinition:
    """Get clock definition by name. Raises KeyError if not found."""
    if name not in CLOCK_DEFINITIONS:
        available = ", ".join(sorted(CLOCK_DEFINITIONS.keys()))
        raise KeyError(f"Clock '{name}' not found. Available: {available}")
    return CLOCK_DEFINITIONS[name]


def search_clocks(
    tissue: str | None = None,
    species: str | None = None,
    output: str | None = None,
) -> dict[str, ClockDefinition]:
    """Filter clock definitions by tissue, species, or output type."""
    results: dict[str, ClockDefinition] = {}
    for name, defn in CLOCK_DEFINITIONS.items():
        if tissue and tissue.lower() not in defn.tissue.lower():
            continue
        if species and species.lower() != defn.species.lower():
            continue
        if output and output.lower() not in defn.output.lower():
            continue
        results[name] = defn
    return results
