begin;
-- recode 1990s
alter table raw90s add column sex int; 
alter table raw90s add column race int;
alter table raw90s add column age5 int;
update raw90s set sex = 
	case 
		when (racesexwtf in (1,3,5,7)) then 1
		when (racesexwtf in (2,4,6,8)) then 2
	end;
update raw90s set race = 
	case 
		when (racesexwtf in (1,2)) then 1
		when (racesexwtf in (3,4)) then 2
		when (racesexwtf in (5,6)) then 3
		when (racesexwtf in (7,8)) then 4
	end;
update raw90s set age5 = 
	case 
		when (agegroup in (0,1)) then 1
		else agegroup
	end;
select * from raw90s where sex is null;
select * from raw90s where race is null;
select * from raw90s where age5 is null;
create table raw90s_sum as 
	select year, fips, sex, age5, sum(pop) as pop
		from raw90s
		group by year, fips, sex, age5;


commit;
