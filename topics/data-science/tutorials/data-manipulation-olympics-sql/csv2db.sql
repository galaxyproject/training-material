.mode csv
.import olympics.tsv.csv olympics
.import olympics_2022.tsv.csv olympics_2022
.import country-information.tsv.csv countries

.mode column
.headers on
.separator ROW "\n"
.nullvalue NULL

update olympics set medal = NULL where medal == "NA";
update olympics set height = NULL where height == "NA";
update olympics set weight = NULL where weight == "NA";
update olympics_2022 set medal = NULL where medal == "NA";
update olympics_2022 set height = NULL where height == "NA";
update olympics_2022 set weight = NULL where weight == "NA";
