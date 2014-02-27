alter table raw00s add column fips text;
update raw00s set fips = to_char(state, 'FM00') || to_char(county, 'FM000');
commit;
