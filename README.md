# Herring spawning migration and bioenergetics 2016 year-class

A simple migration and bioenergetics model for herring spawning migration.

This R-script will replicate the analyses performed in the manuscript.

# Getting the code

```bash
git clone https://github.com/eamousing/herringspawn.git
```

# Software dependencies

The following software versions have been used. Older/newer versions may work but have not been tested.

- R version 4.3.1 (2023-06-16)
- tidyverse version 2.0.0
- ggplot2 version 3.4.3
- patchwork version 1.1.3
- maps version 3.4.1
- geosphere version 1.5.18
- xlsx version 0.6.5
- R.matlab version 3.7.0

# Data dependencies

# Running the analysis

To run the analysis from the terminal on UNIX/Linux:

```bash
Rscript herringspawn_main.R
```

# Analysis output

The results of the migration/bioenergetics model will be placed in `./results` as csv files. In addition, a plot will be produced (`migration.tiff`)summarizing the main results.

# Licence disclaimer

The source code in this project is licensed under the terms of the MIT license. Data provided in the `./data` is **not** subject to the MIT license.
