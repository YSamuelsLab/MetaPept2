perl -ne '@f=split(/,/);$f[2]=~s/ +\+.*//;if($f[2]!~m/^Description/ && $f[1] eq "Linear peptide"){print(">$f[2]\n$f[2]\n")}' epitope_table_export_1668681696.csv > epitope_table_export_1668681696.linear_peptides.fa

perl -ne '@f=split(/\t/);$f[2]=~s/ +\+.*//;if($f[2]!~m/^Description/ && $f[1] eq "Linear peptide"){print(">$f[2]\n$f[2]\n")}' epitope_table_export_1695119936.mice.csv > epitope_table_export_1695119936.mice.fa

