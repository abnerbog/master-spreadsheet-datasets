rule all:
	input: 
		"3_munge/out/master_spreadsheet_formatted.xlsx"

rule fetch_report:
	params:
		api_endpoint = "https://dsp-reports-jbzfw6l52q-uc.a.run.app/csv"
	output:
		out_filename = "1_fetch/out/tmp/discover_report.pkl"
	script:
		"1_fetch/src/report_retrieval.py"

rule clean_report:
	input:
		in_filename = "1_fetch/out/tmp/discover_report.pkl"
	output:
		out_filename = "2_clean/out/tmp/discover_report_cleaned.pkl"
	script:
		"2_clean/src/report_cleaning.py"

rule create_spreadsheet:
	input:
		in_filename = "2_clean/out/tmp/discover_report_cleaned.pkl"
	output:
		out_filename = "3_munge/out/tmp/master_spreadsheet_unformatted.xlsx"
	script:
		"3_munge/src/creating_spreadsheet.py"

rule format_spreadsheet:
	input:
		in_filename = "3_munge/out/tmp/master_spreadsheet_unformatted.xlsx"
	output:
		out_filename = "3_munge/out/master_spreadsheet_formatted.xlsx"
	script:
		"3_munge/src/formatting_spreadsheet.py"
