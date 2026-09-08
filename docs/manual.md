## Running ceif
ceif is a command-line program controlled by arguments. The basic syntax is:

    ceif [OPTION]...

Input data is assumed to be comma-separated values. A different separator can be specified with option -f.
### Options

| Option | Purpose  |
|:----|----|
| -h&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;| Display help and exit|
| -d&nbsp;INTEGER | Number of decimals when printing and saving dimension values. Default is 6|
| -V | Display version and exit|
| -I&nbsp;LIST | LIST is a comma-separated list of field numbers (first = 1) that should be ignored when reading the input file. Ranges can be specified with a dash (e.g., 2-9). The default is to read all fields|
| -U&nbsp;LIST | LIST is a comma-separated list of field numbers (first = 1) that should be processed when reading the input file. Ranges can be specified with a dash (e.g., 2-9). This overrides overlapping values from option -I|
| -X&nbsp;LIST | LIST is a comma-separated list of field numbers (first = 1) that should be treated as text fields when reading the input file. Ranges can be specified with a dash (e.g., 2-9). A hash value in the range 0–32770 is generated from the string value. Note that this is not collision-free and should be used mainly for simple classifications like "yes/no" or "Male/Female/Unknown". Note also that "Yes" and "yes" produce different values|
| -G&nbsp;LIST | LIST is a comma-separated list of dimension attribute indices. A combination of these dimension attributes must yield an outlier score alongside the total score. Ranges can be specified using a dash. Note that these are dimension attribute indices, not input line field indices (first = 1)|
| -t&nbsp;INTEGER&nbsp;&nbsp; | Number of trees to use. Default is 100|
| -s&nbsp;INTEGER | Number of samples for each tree. Default is 256|
| -f&nbsp;CHAR | Field separator for input files|
| -l&nbsp;FILE | File to be used for algorithm training| 
| -a&nbsp;FILE | File to analyze|
| -c&nbsp;FILE | File to categorize|
| -p&nbsp;STRING | Printf-style format to print anomaly data or categorized data. See printing directives below|
| -o&nbsp;FILE | Print output to FILE. Default is stdout|
| -r&nbsp;FILE | Read forest data from file. The file should have been written earlier with option -w|
| -w&nbsp;FILE | Write forest data to FILE. Typically the result of analyzing data from a file with option -l. Data can later be read with option -r|
| -z&nbsp;FILE | Read and write forest data from/to the same file. Forest data is read from FILE and, after processing, written back to FILE|
| -O&nbsp;FLOAT| Outlier score threshold for anomaly detection. Data with a higher or equal score is considered an anomaly and printed using the format given by option -p. Values range from 0.0 to 1.0|
| -O&nbsp;FLOATs| Scaled outlier score for anomaly detection. The analyzed score is scaled to the range 0.0–1.0 using the forest min/max scores. This ensures the best inlier receives a score near 0.0 and the farthest outlier receives a score near 1.0. When given with the categorize option (-c), results are filtered by this threshold: only category results with scores lower than this value are accepted. Values range from 0.0s to 1.0s|
| -O&nbsp;FLOAT%| Outlier score threshold calculated from the sample score distribution, selecting the score value under which FLOAT percent of samples have a lower score. Values range from 0% to 100%|
| -C&nbsp;LIST | List of field numbers to be used as category fields. Default is not to use category fields. Field values are concatenated with colons to form a category string|
| -L&nbsp;LIST | List of field numbers to be used as label fields. Default is not to use label fields. Field values are concatenated with colons to form a label string|
| -F&nbsp;REGEXP | Filter categories using a regular expression. Forests with category strings matching REGEXP are excluded from analysis or categorization. Several options can be specified. If REGEXP is preceded by "-v ", matching is inverted|
| -H | Input data contains a header line that should be ignored. Default is to read all lines|
| -S | Set locale from the environment. Default is the "C" locale|
| -T&nbsp;FLOAT| Generate test data. Test data is generated using the sample set min/max values and test data point interval given by option -i (default is 256). The test data range can be scaled by FLOAT (e.g., 1.0 doubles the test data range). After test data is generated, up to 10,240 sample data points are printed with score 0|
| -i&nbsp;INTEGER| Test data point interval. Larger values produce a denser test dataset|
| -u&nbsp;INTEGER| Accept only unique samples when sampling input data. INTEGER is a value between 0 and 100 (default is 10), representing the percentage of input data rows checked for uniqueness. Value 100 requires all accepted sample data to be unique|
| -m&nbsp;STRING| Printf format for printing floating-point values for samples and sample averages. Default is "%.*f"|
| -j&nbsp;STRING| Printf format for printing all single dimension attribute metrics together. Metrics are formatted using printf directive %m|
| -e&nbsp;CHAR| Value separator when printing sample, sample average, and analyzed data values. Default is comma|
| -M&nbsp;STRING | Print category values, average values, or last update times of forests that were not used in analysis. Optional printf format STRING is used for printing|
| -D&nbsp;INTEGER | Before saving forest data to a file, delete forests that have not been updated within the last INTEGER seconds. If INTEGER is followed by a letter from {Y, M, D, m}, INTEGER is interpreted as years, months, days, or minutes|
| -N&nbsp;STRING | Print input values that are not associated with any category. This can be used for identifying "new" category values. Optional printf format STRING is used for printing|
| -A | Instead of taking samples as individual rows, aggregate new sample values for each forest. Only one new aggregated sample per forest is added for each invocation of option -l|
| -q | Print forest information in human-readable form and exit|
| -y | Print an ASCII density map of forest information and exit|
| -yy | Print an ASCII density map of forest information using a common sample scale for all forests and exit|
| -E | Print samples with their sample scores and exit|
| -k | Remove the sample having the maximum sample score for each non-filtered forest. If specified multiple times, multiple outliers are removed. The updated sample set can be saved with option -w|
| -g&nbsp;FILE | Use FILE as the rc-file instead of ~/.ceifrc. Note that options in FILE override options specified before -g|
| -P | Print a list of correlation coefficients with regression line slopes and y-intercepts for every dimension attribute pair and exit. The correlation coefficient ranges from -1.0 to 1.0|
| -R&nbsp;STRING | Remove all samples for the forest whose category string matches STRING|
| -v&nbsp;STRING | Print average score and other summary statistics calculated from analyzed data using format STRING|
| -Q&nbsp;STRING | Replace input data values using an expression in STRING. STRING is added to the list of expressions. If STRING starts with a hyphen, the expression is removed from the list|


If FILE is "-" then standard input or output is read or written.

Default file format for options -r,-w and -z is JSON. If JSON is not available then CSV format is used. Ceif tries to obey the number of decimals (option -d) when saving data.
If no double formatting support is available, the number of decimals saved is the json library default.

#### Printing directives

| Directive | Meaning |
|----|----|
| %r | Current input file row number|
| %s | Anomaly score|
| %g | Anomaly score for dimensions given by option -G|
| %S | Average anomaly score for analyzed data|
| %n | Number of rows for a forest|
| %o | Number of analyzed rows for a forest. This might be lower than %n if data sampling is enabled (see ANALYZE\_SAMPLING)|
| %h | Number of analyzed rows having a score greater than the outlier score threshold|
| %c | Category string from input data. The original category when categorizing data|
| %C | Forest category string. The best-matching category when categorizing data|
| %l | Label values|
| %d | Separated list of dimension values. Text-based dimensions are printed as text|
| %u | Separated list of dimension values. Text-based dimensions are printed as float|
| %a | Separated list of dimension average values|
| %e | Separated list of dimension attribute scores: ceif analyzes how each attribute affects the total score and assigns each an individual score|
| %m | Separated list of dimension metrics printed attribute by attribute. The printf format for attribute metrics is specified by option -j or the rc-file variable PRINT\_DIMENSION. See rc-file section for printing directives|
| %i | Dimension attribute index (first = 1)|
| %v | Current input row values|
| %x | Outlier score formatted as an RGB hex value (e.g., 127F77). The default color gradient for scores from 0.0 to 1.0 ranges from yellow to red. Colors can be customized in the rc-file. Note that a score of zero is printed as black|
| %X | Outlier score formatted as an RGB hex value for dimensions specified by option -G|
| %t | Timestamp when category was last updated, in human-readable form using current locale|
| %: | Category value separator|
| %. | Label value separator|
| %% | Percent sign|

The value separator for d, u, a, m, e, and v can be specified using option -e.
The category value separator is a semicolon and the label value separator is a dash. These can be customized in the rc-file.

### User rc-file
Default settings can be loaded from the user-specific rc-file `~/.ceifrc`. The file contains variable-value pairs separated by whitespace. Comments start with `#`. 
These values can be overridden by command options and loaded forest data (read via option -r). Note that a file specified with option -g overrides options specified prior to -g.

The following variables are supported:

| Variable | Meaning | Default Value |
|----|----|----|
|SAMPLES|Number of samples taken for each tree; same effect as option -s|256|
|TREES|Number of trees for each forest; same effect as option -t|100|
|DECIMALS|Number of decimals used when saving forest data; also affects printing of sample values (option -d)|6|
|AUTO\_SCALE|Scale sample values before analyzing the forest: 1 = yes, 0 = no|1|
|CATEGORY\_SEPARATOR|Character used as a separator when concatenating category fields|;|
|LABEL\_SEPARATOR|Character used as a separator when concatenating label fields|-|
|OUTLIER\_SCORE|Outlier score threshold for analysis; accepts the same formats as option -O ("max", "average", float 0..1, or scaled float 0s..1s)|0.65|
|MAX\_SAMPLES|Maximum number of samples for each forest|Calculated as number\_of\_trees * number\_of\_samples\_per\_tree|
|NEAREST|Score is adjusted by distance to nearest sample point in leaf nodes: 1 = yes, 0 = no|1|
|ANALYZE\_SAMPLING|If analyzed data is impractically large, stream sampling can be enabled. When the analyzed row count reaches this threshold, sampling begins using reservoir sampling. The number of analyzed rows is estimated as k * (ln(x/k) + 1), where k is this parameter and x is total rows|0 (disabled)|
|DEBUG|Print debug messages: 1 = yes, 0 = no|0|
|PRINT\_DIMENSION|Printf format string for directive %m. May contain directives %d, %a, %e, and %i|""|
|DIM\_PRINT\_WIDTH|Attribute metric column width when printing forest info with option -q|25|
|CLUSTER\_SIZE|ceif identifies data clusters by selecting samples with the lowest scores and counting adjacent samples. Cluster radius is calculated by finding the distance from the lowest-scoring sample to the most distant sample, multiplied by this parameter (range: 0 to 1)|0.125|
|LOW_RGB_COLOR|RGB color code for score 0 (%x directive). Given as hex string (e.g., 0xffff00)|0xffff00 (yellow)|
|HIGH_RGB_COLOR|RGB color code for score 1 (%x directive)|0xff0000 (red)|

Example rc-file:

    # My default values
    TREES 200
    CATEGORY_SEPARATOR +
    LABEL_SEPARATOR   .
    OUTLIER_SCORE 0.65s
    PRINT_DIMENSION "<value>%d</value>\n<score>%e</score>\n"

    # End of file

#### Deprecated Parameters

| Variable | Status / Description |
|:---|:---|
| `CENTROID_THRESHOLD`<br>*(or `CENTROID_TRESSHOLD`)* | **Deprecated.** In earlier versions, this specified the tree depth ratio threshold at which node splitting switched from random perturbations to sample centroids for tree balancing. Tree construction now uses pairwise local sample interpolation (SCiForest / Pairwise EIF) across all depths, which naturally stays bounded within local cluster geometry and cleanly separates samples without requiring centroid balancing. This parameter and the associated API (`set_centroid_tresshold`) are preserved for backwards compatibility but silently ignored. |
| `AUTO_WEIGTH` | **Deprecated.** Legacy spelling and alias for `AUTO_SCALE`. |

### Examples

#### Learn and analyze the same file
Learn and analyze `data.csv` using fields 1-3 (ignoring fields 4-100), and print anomalies having a score of 0.6 or greater:

    ceif -l data.csv -a data.csv -I4-100 -O0.6

#### Learn and write forest data to data.f with scaled score threshold 0.5s:

    ceif -l data.csv -w data.f -I4-100 -O0.5s

#### Read forest data from data.f and analyze data.csv
Note that parameters such as -O are saved inside the forest file and do not need to be specified again:

    ceif -r data.f -a data.csv

#### Read forest data from data.f and add more samples from data2.csv 
Note that -w must be specified in order to save the updated forest model:

    ceif -r data.f -l data2.csv -w data.f

Option -z reads and writes back to the same file:

    ceif -z data.f -l data2.csv 

#### Learn and write categorized forest data 
Use field number 5 as the category key:

    ceif -l data.csv -w data.f -I4-100 -O0.6 -C 5

#### Categorize data from data.csv using the forest data from the previous example
Print the analyzed category match (%C) and the dimension values (%d) for each input row from `data.csv`:

    ceif -r data.f -c data.csv -p "%C %d"

#### Generate a test dataset using a forest data file
Generate a synthetic evaluation dataset around sample data points by enlarging the range by a factor of 1. Each dimension range evaluates 512 test points. Test data values and outlier scores formatted as RGB hex values are separated by semicolons and written to `plot_data.csv`:

    ceif -r data.f -T1 -i512 -e";" -p"%d;0x%x" -o plot_data.csv

#### Data value aggregation
If daily or hourly summaries are to be analyzed, `ceif` can be invoked with the -A option, where each forest represents an aggregation key. 
Example CSV file of hourly internet traffic samples:

    $ cat traffic.csv
    Hour,Type,Inbytes,Outbytes
    12,In,123,1444
    12,Out,423,1644
    13,In,123,1444
    13,In,823,44
    14,In,13,1414
    14,In,9123,1444
    14,In,123,1443
    15,Out,12423,1644
    16,Out,423,16
    16,Out,493,1044
    17,Out,433,1644
    18,Out,493,1644

Generating hourly and traffic-direction summaries, where the first two fields serve as category keys (-C1-2):

    ceif -l traffic.csv -C1-2 -H -A -d0 -w traffic.ceif

Forest file contents below; for each hour and traffic type, byte counts are summed into aggregated samples:

    $ cat traffic.ceif
    G;2;"";"%s %v";100;256;"1-2";",";1;0.750000;10.000000;"";"";8;"";0;0;"";",";1;1
    F;"12:In";0.000000;0;1;1582453164
    S;123|1444
    F;"12:Out";0.000000;0;1;1582453164
    S;423|1644
    F;"13:In";0.000000;0;1;1582453164
    S;946|1488
    F;"14:In";0.000000;0;1;1582453164
    S;9259|4301
    F;"15:Out";0.000000;0;1;1582453164
    S;12423|1644
    F;"16:Out";0.000000;0;1;1582453164
    S;916|1060
    F;"17:Out";0.000000;0;1;1582453164
    S;433|1644
    F;"18:Out";0.000000;0;1;1582453164
    S;493|1644

When this is repeated on subsequent days, a new aggregated sample row is added to each forest.
Once sufficient historical baseline samples have been gathered, daily traffic can be evaluated for anomalies using option -a.

#### Using input data expressions (-Q)
Input data fields can be dynamically transformed before processing. Expressions support standard `tinyexpr` syntax (see https://github.com/codeplea/tinyexpr). Input fields are referenced using the `$n` notation, where `n` is the 1-based index of the input field.

The expression format is:
```
$n=<expr>
```
or with decimal rounding:
```
$n=<expr>:d
```
Where:
- `$n` is the field to be modified.
- `<expr>` is any valid `tinyexpr` expression that may include `$n` references to input fields.
- `:d` is the optional number of decimal digits to retain.

Examples:

| Case | Expression |
|:---|:---|
| Divide the 4th field by 13, retaining two decimals | `'$4=$4/13:2'` |
| Multiply the first field by the second field | `'$1=$1*$2'` |
| Divide the second field by 60 and floor to an integer | `'$2=floor($2/60)'` |
| Assign a constant value to field 10 | `'$10=3486'` |

All active expressions are persisted inside the saved forest metadata. To remove a previously saved expression, prefix it with a hyphen (e.g., `-Q '-$4=$4/10'`) and re-save the forest (`-w` or `-z`).

Example command line:
```bash
ceif -Q '$4=$4/13:2' -Q '$2=floor($2/60)' -Q '-$4=$4/10' ...
```
