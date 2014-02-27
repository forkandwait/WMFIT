W. Webb Sprague
webb.sprague@gmail.com
2012-06-01
Copyright blah, make sure you reference me.

--- Overview ---

The files in this directory are age and sex population forecasts for most of
counties in the USA, generated using the method described in the article
"Automatic parametrization of age/ sex Leslie matrices for human populations"
(http://arxiv.org/abs/1203.2313)


--- File Descriptions ---

There are three data files in the zip package:  

1.  "allcntysrpt-TIMESTAMP-fcst.txt" 

    This file holds population counts for all the counties for which data was
    available and which converged when we tried to fit a Leslie transition
    matrix using Wood's Method.  It also holds cross-validation forecasts.

    The columns are named at the top of the file according the year and whether
    they are empirical ("E") or forecasted ("F").  F2010 is forecasted from
    E2000 as the jump off year, using 1980 to 2000 age sex population counts to
    generate the transition matrix. F2015 through F2040 are forecasted from 2010
    using 1980 through 2010 input to generate the matrix.

    Measures of forecast error can be calculated by comparing F2010 and E2010,
    (forecasted 2010 counts and empirical 2010 counts, respectively).  Assuming
    a single county is stored in a variable "p", In Octave/ Matlab, to calculate
    MAPE:

	mape = mean (abs(p(:,2) .- p(:,3)) ./ p(:,3)); 


2.  "allcntysrpt-TIMESTAMP-m00.txt" 

    This file holds Leslie matrixes for 2010 cross validation forecasts.  These
    matrices were generated from 1980 to 2000 input, then used to forecast to
    2010; the results were compared to empirical 2010 data to analyze error.

    The user should be able to re-derive F2010 by simple linear algebra and the
    jump off data in the file described above:

    F2010 = (m00 ^ 2) * F2000

    To visualize the difference between forecast and empirical populations in
    Octave/ Matlab, load the matrix into a variable called "m00", load the
    population into a variable called "p", and then run the following:

    plot([m00^2 * p(:,1), p(:,3)]);
    

3.  "allcntysrpt-TIMESTAMP-m10.txt" 

    This file holds Leslie matrices that generate the forecasts for 2015 through
    2040 (F2015, F2020, etc.).  These matrices can be used to re-derive the
    included forecast numbers using the same equation as above, varying the
    exponent depending on how many five year periods to project:

    F2015 = (m10 ^ 1) * E2000

    F2020 = (m10 ^ 2) * E2000

    Etc.

    In Octave/ Matlab:  plot([m10 * p(:,3), p(:,4)]);


--- Notes ---

About missing FIPS: Some FIPS were dropped from the analysis, due to naming
nconsistencies, lack of convergence, other random stuff.

About numerical differences: errors of +/- 15.0 may have crept in from rounding,
numerical leprechauns, etc.  The Leslie matrix is formatted for readability, and
the dropping of decimal can have a large effect.

About structured text files: The text files have been written so that they are
parseable by programs. There is a FIPS, then the data for a county, then an
empty line.  One can use regular expressions to search for the fips
(e.g. "^53001" to find Adams County in Washington).  The three data files are
also "parallel" to each other, so the 2000th record (records are defined by an
empty line separator) in all three files should refer to the same county in all
three files (XXX county...).

About importing data: In Octave/ Matlab, one can cut and paste the matrices by typing
"[" at the terminal, pasting the data, typing a "]", and hitting enter.

About delimiters: All the data files are space delimited.  With LibreOffice
Calc, choose "merge delimiters" in the input dialogue.

About filenames: The "TIMESTAMP" reference above is formatted as YYYYmmDDTHHMMSS
-- YYYY represents a four digit year, mm a two digit month, DD a two digit day,
T is a separator, HH is a two digit hour, MM a two digit minute, and SS a two
digit second.  All of these can have leading zeros, so May would be written as
"05", for example.

About Octave: it rocks, and can be downloaded for free.  Google is your friend.
Then open an Octave terminal, and follow a Matlab tutorial.

--- Closing --- 

I am available for contract work, or will do whatever if it is interesting enough.

Send corrections and such to me.
