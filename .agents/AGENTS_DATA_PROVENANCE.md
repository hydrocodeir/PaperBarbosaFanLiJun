# AGENTS_DATA_PROVENANCE.md

## Role

Maintain complete scientific provenance for every dataset used or generated in the project.

The objective is to make it possible to answer:

> Where did this value come from, how was it transformed, and which scientific decision produced it?

---

## Scope

Track:

- raw data;
- downloaded data;
- user-provided data;
- derived data;
- intermediate data;
- model inputs;
- model outputs;
- calibration data;
- validation data;
- masks;
- boundaries;
- metadata;
- lookup tables;
- manually entered values.

---

## Core Principle

No scientifically important dataset should appear in the analytical workflow without documented origin.

---

## Dataset Record

For each dataset document:

```markdown
## Dataset: <name>

### Identity
- Official name:
- Short name:
- Version:
- Provider:
- DOI:
- Landing page:
- License:

### Acquisition
- Acquisition date:
- Acquisition method:
- API / download / user supplied / generated:
- Original format:

### Coverage
- Spatial extent:
- Spatial resolution:
- Temporal extent:
- Temporal resolution:

### Variables
| Variable | Meaning | Unit | Original Unit | Notes |
|---|---|---|---|---|

### Spatial Metadata
- CRS:
- EPSG:
- Grid type:
- Longitude convention:
- Latitude orientation:

### Missing Data
- Missing-value code:
- Missingness treatment:

### Quality Control
- QC flags:
- Filters:
- Known limitations:

### Preprocessing
1. ...
2. ...

### Derived Products
- ...

### Validation
- ...

### Provenance Status
- Verified / Partially Verified / Unverified
```

---

## Raw Data Preservation

Whenever practical:

- preserve raw data unchanged;
- do not overwrite originals;
- separate raw and processed products;
- preserve original metadata;
- preserve checksum or file fingerprint when useful.

Never silently modify source data.

---

## Transformation Chain

Every derived dataset must document its parent dataset and transformations.

Example:

```text
ERA5 hourly temperature
    ↓ unit conversion
Kelvin → Celsius
    ↓ daily aggregation
Daily Tmax
    ↓ spatial clipping
Study-area Tmax
    ↓ index computation
TXx
```

The lineage must remain traceable.

---

## Unit Tracking

For every variable record:

- source unit;
- processing unit;
- reporting unit;
- conversion formula.

Never infer unit solely from numerical magnitude when metadata are available.

---

## Spatial Provenance

For geospatial data document:

- CRS;
- datum;
- grid;
- resampling;
- reprojection;
- clipping;
- interpolation;
- rasterization;
- vectorization;
- cell alignment.

When combining grids, identify which dataset defined the target grid.

---

## Temporal Provenance

Document:

- timezone;
- calendar;
- leap-day treatment;
- aggregation;
- resampling;
- temporal alignment;
- missing periods;
- climatological baseline.

For climate data, record non-Gregorian calendars where relevant.

---

## Data Versioning

If a dataset is updated:

- retain old provenance;
- record new version;
- describe differences;
- identify analyses affected.

Do not silently replace one release with another.

---

## Derived Variables

For each derived variable record:

- equation;
- source variables;
- units;
- thresholds;
- window;
- aggregation;
- implementation reference.

---

## External Data Verification

Before treating external data as authoritative, verify when possible:

- provider;
- dataset name;
- version;
- variable definition;
- resolution;
- period;
- DOI or official documentation.

Never fabricate provenance fields.

---

## Manual Data Entry

If values are manually entered, record:

- source;
- person or document basis if appropriate;
- transcription method;
- verification;
- date entered.

Manual values must not masquerade as automatically generated data.

---

## Data Quality Flags

Classify concerns:

- Missing
- Suspect
- Corrected
- Imputed
- Interpolated
- Excluded
- Unverified

---

## Provenance Table

Maintain a master table:

```markdown
| Dataset | Version | Source | Role | Raw/Derived | Verified | Key Transformations |
|---|---|---|---|---|---|---|
```

---

## Scientific Output Traceability

Every important result should be traceable to:

1. output value;
2. output table/figure;
3. analysis product;
4. processing step;
5. source dataset.

---

## Non-Negotiable Rules

1. Never overwrite raw data without explicit instruction.
2. Never invent missing metadata.
3. Never assume CRS, units, calendar, or resolution when they can be verified.
4. Never lose the transformation chain.
5. Never present derived data as raw observations.
6. Never use undocumented data in final scientific conclusions.
