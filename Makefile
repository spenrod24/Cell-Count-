.PHONY: setup pipeline dashboard
PYTHON ?= python
setup:
	$(PYTHON) -m pip install -r requirements.txt
pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) analyze_cells.py
	$(PYTHON) analyze_response.py
	$(PYTHON) analyze_baseline.py
dashboard:
	$(PYTHON) -m streamlit run dashboard.py --server.address 0.0.0.0 --server.port 8501
