-- Not all fips codes are in all files ...  TODO deal with changes
-- (Dade County!) instead of clobbering.  Also, boundary changes might
-- happen even with constant names.

begin;
drop schema fips if exists;
create schema fips;

create table fips.f70s as select distinct fips from raw70s;
create table fips.f80s as select distinct fips from raw80s;
create table fips.f90s as select distinct fips from raw90s;
create table fips.f00s as select distinct fips from raw00s;

create table fips.fips_union as (
	select fips from fips.f70s
		union
	select fips from fips.f80s
		union
	select fips from fips.f90s
		union
	select fips from fips.f00s
);


create table fips.fips_inter as (
	select fips from fips.f70s
		intersect
	select fips from fips.f80s
		intersect
	select fips from fips.f90s
		intersect
	select fips from fips.f00s
);

create table fips.wtf as (
	select fips from fips.fips_union 
		except
	select fips from fips.fips_inter);

delete from raw70s_sum a where exists (select 1 from fips.wtf b where a.fips=b.fips);
delete from raw80s_sum a where exists (select 1 from fips.wtf b where a.fips=b.fips);
delete from raw90s_sum a where exists (select 1 from fips.wtf b where a.fips=b.fips);
delete from raw00s a where exists (select 1 from fips.wtf b where a.fips=b.fips);

select count(*) from (select distinct fips from raw70s_sum) a;
select count(*) from (select distinct fips from raw80s_sum) a;
select count(*) from (select distinct fips from raw90s_sum) a;
select count(*) from (select distinct fips from raw00s) a;

commit;

/*
webbs=# select * from fips.wtf order by fips;  
 fips
-------
 02010
 02013
 02016
 02030
 02040
 02068
 02080
 02105
 02120
 02122
 02140
 02160
 02164
 02185
 02188
 02190
 02195
 02198
 02200
 02201
 02210
 02230
 02231
 02232
 02250
 02260
 02261
 02275
 02280
 02282
 04012
 08014
 12025
 12086
 15005
 29186
 29193
 30113
 35006
 46131
 51123
 51560
 51683
 51685
 51735
 51780
(46 rows)

*/

