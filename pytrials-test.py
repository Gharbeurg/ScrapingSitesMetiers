from pytrials.client import ClinicalTrials
import csv
import os

ct = ClinicalTrials()

# Get the NCTId, Condition and Brief title fields from 500 studies related to Coronavirus and Covid, in csv format.
corona_fields = ct.get_study_fields(
    search_expr="AREA[FunderTypeSearch]industry AND AREA[StartDate]RANGE[08/01/2022, 01/01/2023] AND AREA[Phase]Phase 3",
    fields=["NCTId", "ResponsiblePartyInvestigatorAffiliation", "Condition", "ConditionMeshTerm", "BriefTitle", "OfficialTitle", "BriefSummary", "DetailedDescription", "Keyword", "EnrollmentCount", "PrimaryOutcomeMeasure", "LeadSponsorName", "Phase", "StartDate", "CompletionDate", "LocationCountry"],
    max_studies=1000,
    fmt="csv",
)
if os.path.exists('C:/PYTHON/.data/trials.csv'):
    os.remove('C:/PYTHON/.data/trials.csv')

with open('C:/PYTHON/.data/trials.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    writer.writerows(corona_fields)