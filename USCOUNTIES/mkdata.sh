set -eu

## Get 70s

SQL="
    select age0
           , age5
           , age10
           , age15
           , age20
           , age25
           , age30
           , age35
           , age40
           , age45
           , age50
           , age55
           , age60
           , age65
           , age70
           , age75
           , age80
           , age85
        from raw70s_sum
        where yr=1970
        order by fips, sex;
"
psql -qAt --port=5433 --dbname=intercensal --command="$SQL" | gawk -f w2l.awk > d70.txt

SQL="
    select age0
           , age5
           , age10
           , age15
           , age20
           , age25
           , age30
           , age35
           , age40
           , age45
           , age50
           , age55
           , age60
           , age65
           , age70
           , age75
           , age80
           , age85
        from raw70s_sum
        where yr=1975
        order by fips, sex;
"
psql -qAt --port=5433 --dbname=intercensal --command="$SQL" | gawk -f w2l.awk > d75.txt

## Get 80s

SQL="
    select age0
           , age5
           , age10
           , age15
           , age20
           , age25
           , age30
           , age35
           , age40
           , age45
           , age50
           , age55
           , age60
           , age65
           , age70
           , age75
           , age80
           , age85
        from raw80s_sum
        where yr=1980
        order by fips, sex;
"
psql -qAt --port=5433 --dbname=intercensal --command="$SQL" | awk -f w2l.awk > d80.txt

SQL="
    select age0
           , age5
           , age10
           , age15
           , age20
           , age25
           , age30
           , age35
           , age40
           , age45
           , age50
           , age55
           , age60
           , age65
           , age70
           , age75
           , age80
           , age85
        from raw80s_sum
        where yr=1985
        order by fips, sex;
"
psql -qAt --port=5433 --dbname=intercensal --command="$SQL" | awk -f w2l.awk > d85.txt

## Get 90s

psql -qAt --port=5433 --dbname=intercensal --command="select pop from raw90s_sum where yr=90 order by fips, sex, age5;" > d90.txt
psql -qAt --port=5433 --dbname=intercensal --command="select pop from raw90s_sum where yr=95 order by fips, sex, age5;" > d95.txt

## Get 2000s

SQL="
    select popestimate2000 from raw00s 
        where sex in (1,2) and agegrp <> 0 
        order by fips, sex, agegrp;
"
psql -qAt --port=5433 --dbname=intercensal --command="$SQL" > d00.txt

SQL="
    select popestimate2005 from raw00s 
        where sex in (1,2) and agegrp <> 0 
        order by fips, sex, agegrp;
"
psql -qAt --port=5433 --dbname=intercensal --command="$SQL" > d05.txt

SQL="
    select popestimate2010 from raw00s 
        where sex in (1,2) and agegrp <> 0 
        order by fips, sex, agegrp;
"
psql -qAt --port=5433 --dbname=intercensal --command="$SQL" > d10.txt

## Get FIPS/sex/age thing to match
SQL="select fips, sex, agegrp from raw00s where sex in (1,2) and agegrp <> 0 order by fips, sex, agegrp;";
psql -qAt --field-separator="	" --port=5433 --dbname=intercensal --command="$SQL" > dkeys.txt

## Make a big matrix file
paste dkeys.txt d70.txt d75.txt d80.txt d85.txt d90.txt d95.txt d00.txt d05.txt d10.txt > dall.txt

echo "Check number of lines of everything!"
wc -l d*.txt

# push intermediates into JUNK
mv dkeys.txt d70.txt d75.txt d80.txt d85.txt d90.txt d95.txt d00.txt d05.txt d10.txt JUNK
