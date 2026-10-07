# CheckMol functional group census

Counts occurrences of 204 functional group types in a molecule, spanning alcohols, carbonyls, carboxylic acids, amines, halides, heterocycles and many more. Haider released checkmol as an open-source command-line program alongside the MolDB structure-database package, and it perceives groups through explicit chemical rules rather than a learned model, so every count traces to a defined group and stays interpretable, unlike a hashed fingerprint. Molecules that cannot be parsed return empty values rather than zeros, which separates failure from genuine absence.

This model was incorporated on 2026-06-02.Last packaged on 2026-06-02.

## Information
### Identifiers
- **Ersilia Identifier:** `eos5f0j`
- **Slug:** `checkmol-functional-groups`

### Domain
- **Task:** `Representation`
- **Subtask:** `Featurization`
- **Biomedical Area:** `Any`
- **Target Organism:** `Any`
- **Tags:** `Fingerprint`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `204`
- **Output Consistency:** `Fixed`
- **Interpretation:** Counts of each of 204 functional group types detected in the molecule.

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| cation | integer | high | Number of cation occurrences in the molecule |
| anion | integer | high | Number of anion occurrences in the molecule |
| carbonyl_compound | integer | high | Number of carbonyl compound occurrences in the molecule |
| aldehyde | integer | high | Number of aldehyde occurrences in the molecule |
| ketone | integer | high | Number of ketone occurrences in the molecule |
| thiocarbonyl_compound | integer | high | Number of thiocarbonyl compound occurrences in the molecule |
| thioaldehyde | integer | high | Number of thioaldehyde occurrences in the molecule |
| thioketone | integer | high | Number of thioketone occurrences in the molecule |
| imine | integer | high | Number of imine occurrences in the molecule |
| hydrazone | integer | high | Number of hydrazone occurrences in the molecule |

_10 of 204 columns are shown_
### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos5f0j](https://hub.docker.com/r/ersiliaos/eos5f0j)
- **Docker Architecture:** `AMD64`, `ARM64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos5f0j.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos5f0j.zip)

### Resource Consumption
- **Model Size (Mb):** `4`
- **Environment Size (Mb):** `518`
- **Image Size (Mb):** `540.69`

**Computational Performance (seconds):**
- 10 inputs: `34.62`
- 100 inputs: `24.9`
- 10000 inputs: `58.81`

### References
- **Source Code**: [https://homepage.univie.ac.at/norbert.haider/cheminf/cmmm.html](https://homepage.univie.ac.at/norbert.haider/cheminf/cmmm.html)
- **Publication**: [https://doi.org/10.3390/molecules15085079](https://doi.org/10.3390/molecules15085079)
- **Publication Type:** `Peer reviewed`
- **Publication Year:** `2010`
- **Ersilia Contributor:** [miquelduranfrigola](https://github.com/miquelduranfrigola)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [GPL-3.0-or-later](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos5f0j
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos5f0j
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
