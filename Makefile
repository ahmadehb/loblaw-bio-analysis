# Makefile — graders run 'make setup', 'make pipeline', 'make dashboard'.

PYTHON ?= python3
PIP    ?= pip

.PHONY: setup pipeline dashboard clean

# Install all dependencies
setup:
	$(PIP) install --upgrade pip
	$(PIP) install -r requirements.txt

# Run full data pipeline end-to-end:
#   1. load_data.py -> initialize SQLite DB and load CSV
#   2. src/analysis.py -> per-sample frequency table
#   3. src/statistics.py -> responder vs non-responder stats
#   4. src/subset_analysis.py -> baseline subset queries
pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) src/analysis.py
	$(PYTHON) src/statistics.py
	$(PYTHON) src/subset_analysis.py

# Launch Streamlit dashboard on the default port
# In GitHub Codespaces, the port is auto-forwarded to a preview URL
dashboard:
	streamlit run dashboard/app.py --server.headless true --server.address 0.0.0.0

# Remove generated artifacts (DB and outputs/) so the pipeline can re-run cleanly
clean:
	rm -f loblaw.db
	rm -rf outputs/
