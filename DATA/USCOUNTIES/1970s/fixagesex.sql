begin;

-- recode 1970s and make sum table
alter table raw70s add column sex int; 
alter table raw70s add column race int;
update raw70s set sex = 
	case 
		when (racesexwtf ~* ' male') then 1
		when (racesexwtf ~* ' female') then 2
	end;
update raw70s set race = 
	case 
		when (racesexwtf ~* 'white') then 1
		when (racesexwtf ~* 'black') then 2
		when (racesexwtf ~* 'other') then 3
	end;
select * from raw70s where sex is null;
select * from raw70s where race is null;


create table raw70s_sum as 
	select yr, fips, sex
			, sum(age0) as age0
			, sum(age5) as age5
			, sum(age10) as age10
			, sum(age15) as age15
			, sum(age20) as age20
			, sum(age25) as age25
			, sum(age30) as age30
			, sum(age35) as age35
			, sum(age40) as age40
			, sum(age45) as age45
			, sum(age50) as age50
			, sum(age55) as age55
			, sum(age60) as age60
			, sum(age65) as age65
			, sum(age70) as age70
			, sum(age75) as age75
			, sum(age80) as age80
			, sum(age85) as age85
		from raw70s
		group by yr, fips, sex;


commit;
