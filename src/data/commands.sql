-- SQL File to improve query speeds
-- Run these commands if ever the tables get recreated and there is slowdown

-- Create index
create index outcomes_idx on outcomes (Alias, Gene, fraction_per_barcode)   

-- list indexes
pragma index_list('outcomes');

-- explain query plan
explain query plan
select Barcode, Gene, Alias, outcome, 
        fraction_per_barcode, mutEvents, countEvents 
        from outcomes where Alias in ("MB01", "MB02", "MB03", "MB04", "MB05", "MB06") and 
        Gene is not "Empty" and Gene is not "Non-targeting" and fraction_per_barcode > 0.01
