#!/usr/bin/perl

use strict;
use Math::Complex;

require "$ENV{PERL_HOME}/Lib/load_args.pl";
require "$ENV{PERL_HOME}/Lib/libfile.pl";
require "$ENV{PERL_HOME}/Lib/format_number.pl";
require "$ENV{PERL_HOME}/Sequence/sequence_helpers.pl";

if ($ARGV[0] eq "--help")
{
  print STDOUT <DATA>;
  exit;
}

my $file_ref;
my $file = $ARGV[0];
if (length($file) < 1 or $file =~ /^-/) 
{
  $file_ref = \*STDIN;
}
else
{
  open(FILE, $file) or die("Could not open file '$file'.\n");
  $file_ref = \*FILE;
}

my %args = load_args(\@ARGV);

my $beginning = get_arg("b", 0, \%args);
my $column_location = get_arg("c", -1, \%args);
my $start_counter = get_arg("n", 0, \%args);
my $column_string = get_arg("s", "", \%args);
my $empty_string = get_arg("e", 0, \%args);
my $divide_columns_string = get_arg("d", "", \%args);
my $multiply_columns_string = get_arg("m", "", \%args);
my $subtract_columns_string = get_arg("u", "", \%args);
my $add_columns_string = get_arg("a", "", \%args);
my $min_columns_string = get_arg("min", "", \%args);
my $max_columns_string = get_arg("max", "", \%args);
my $count_columns_string = get_arg("count", "", \%args);
my $average_columns_string = get_arg("ave", "", \%args);
my $std_columns_string = get_arg("std", "", \%args);
my $stderr_columns_string = get_arg("se", "", \%args);
my $qauntile_range_string = get_arg("quant", "0,1", \%args);
my $concat_columns_string = get_arg("cat","", \%args);
my $delim = get_arg("delim","", \%args);

my $add_file = get_arg("f", "", \%args);
my $significant_nums = get_arg("sn", "2", \%args);

my $ignore_zero_entries = get_arg("nozero", 0, \%args);

$column_string =~ s/\"//g;

my @divide_columns = split(/\,/, $divide_columns_string);
my @concat_columns = split(/\,/, $concat_columns_string);
my @multiply_columns = &ParseRanges($multiply_columns_string);
my @subtract_columns = split(/\,/, $subtract_columns_string);
my @add_columns = &ParseRanges($add_columns_string);
my @min_columns = &ParseRanges($min_columns_string);
my @max_columns = &ParseRanges($max_columns_string);
my @count_columns = &ParseRanges($count_columns_string);
my @average_columns = &ParseRanges($average_columns_string);
my @std_columns = &ParseRanges($std_columns_string);
my @stderr_columns = &ParseRanges($stderr_columns_string);
my @qauntile_range = split(/\,/, $qauntile_range_string);

my @columns;
if (length($add_file) > 0)
{
  open(FILE2, "<$add_file");
  while(<FILE2>)
  {
    chomp;

    push(@columns, $_);
  }
}

my $row_counter = 0;

my $mytest= length($add_file);

while (<$file_ref>)
{
  chomp;
  my $str = "";
  if (length($add_file) > 0) 
{ 
      $str = $columns[$row_counter];
}
  elsif (length($column_string) > 0) { $str = $column_string; }
  elsif ( $empty_string ) { $str = ""; }
  elsif (length($divide_columns_string) > 0)
  {
    my @row = split(/\t/);

    $str = (length($row[$divide_columns[0]]) > 0 and length($row[$divide_columns[1]]) > 0 and $row[$divide_columns[1]] != 0) ? &format_number($row[$divide_columns[0]] / $row[$divide_columns[1]], $significant_nums) : "";
  }
  elsif (length($concat_columns_string) > 0)
{
    my @row = split(/\t/);
    $str="";
    my $k=0;
    for($k=0; $k < @concat_columns-1; $k++)
    {
	$str = $str . $row[$concat_columns[$k]] . $delim;
    }
    $str = $str . $row[$concat_columns[$k]];
}

  elsif (length($multiply_columns_string) > 0)
  {
    my @row = split(/\t/);
    $str = 1;
    my $not_all_empty = 0;
    for(my $k=0; $k < @multiply_columns; $k++)
    {
       next if $ignore_zero_entries and $row[$multiply_columns[$k]] == 0;

       if (length($row[$multiply_columns[$k]]) > 0)
       {
	  $not_all_empty = 1;
	  $str = $str * $row[$multiply_columns[$k]];
       }
    }
    $str = ($not_all_empty == 1) ? &format_number($str, $significant_nums) : "";
  }
  elsif (length($subtract_columns_string) > 0)
  {
    my @row = split(/\t/);

    $str = (length($row[$subtract_columns[0]]) > 0 and length($row[$subtract_columns[1]]) > 0) ? &format_number($row[$subtract_columns[0]] - $row[$subtract_columns[1]], $significant_nums) : "";
  }
  elsif (length($add_columns_string) > 0)
  {
    my @row = split(/\t/);
    $str = 0;
    my $not_all_empty = 0;
    for(my $k=0; $k < @add_columns; $k++)
    {
       next if $ignore_zero_entries and $row[$add_columns[$k]] == 0;

       if (length($add_columns[$k]) > 0)
       {
	  $not_all_empty = 1;
	  $str = $str + $row[$add_columns[$k]];
       }
    }
    $str = ($not_all_empty == 1) ? &format_number($str, $significant_nums) : "";
  }
  elsif (length($min_columns_string) > 0)
  {
    my @row = split(/\t/);

    $str = 100000000000000000;
    my $not_all_empty = 0;

    for(my $k=0; $k < @min_columns; $k++)
    {
       next if $ignore_zero_entries and $row[$min_columns[$k]] == 0;

       if (length($row[$min_columns[$k]]) > 0)
       {
	  $not_all_empty = 1;
	  $str = ($str < $row[$min_columns[$k]]) ? $str : $row[$min_columns[$k]];
       }
    }
    $str = ($not_all_empty == 1) ? &format_number($str, $significant_nums) : "";
  }
  elsif (length($max_columns_string) > 0)
  {
    my @row = split(/\t/);

    $str = -100000000000000000;
    my $not_all_empty = 0;

    for(my $k=0; $k < @max_columns; $k++)
    {
       next if $ignore_zero_entries and $row[$max_columns[$k]] == 0;

       if (length($row[$max_columns[$k]]) > 0)
       {
	  $not_all_empty = 1;
	  $str = ($str > $row[$max_columns[$k]]) ? $str : $row[$max_columns[$k]];
       }
    }
    $str = ($not_all_empty == 1) ? &format_number($str, $significant_nums) : "";
  }
  elsif (length($count_columns_string) > 0)
  {
    my @row = split(/\t/);

    $str = 0;

    for(my $k=0; $k < @count_columns; $k++)
    {
       next if $ignore_zero_entries and $row[$count_columns[$k]] == 0;

       if (length($row[$count_columns[$k]]) > 0)
       {
	  $str = $str + 1;
       }
    }
  }
  elsif (length($average_columns_string) > 0)
  {
    my @row = split(/\t/);

    $str = "";
    my $not_all_empty = 0;
    my @tmp_list = ();

    for(my $k=0; $k < @average_columns; $k++)
    {
       next if $ignore_zero_entries and $row[$average_columns[$k]] == 0;

       if (length($row[$average_columns[$k]]) > 0)
       {
	  $not_all_empty = 1;
	  push(@tmp_list,$row[$average_columns[$k]]);
       }
    }

    if ($not_all_empty == 1)
    {
       if (($qauntile_range[0] > 0) or ($qauntile_range[1] < 1))
       {
	  @tmp_list = sort {$a <=> $b} @tmp_list;
       }
       my $first = int($qauntile_range[0] * @tmp_list);
       $first = ($first < 0) ? 0 : $first;
       my $last = int($qauntile_range[1] * @tmp_list);
       $last = ($last > (@tmp_list - 1)) ? (@tmp_list - 1) : $last;
       $first = ($first > $last) ? $last : $first;

       if (($last >= $first) and ($last < @tmp_list))
       {
	  $str = 0;
	  for (my $k = $first; $k <= $last; $k++)
	  {
	     $str = $str + $tmp_list[$k];
	  }
	  $str = $str / ($last - $first +1);
       }
    }
    $str = (($not_all_empty == 1) and (length($str) > 0)) ? &format_number($str, $significant_nums) : "";
  }
 elsif (length($std_columns_string) > 0)
  {
    my @row = split(/\t/);

    $str = "";
    my $not_all_empty = 0;
    my @tmp_list = ();

    for(my $k=0; $k < @std_columns; $k++)
    {
       next if $ignore_zero_entries and $row[$std_columns[$k]] == 0;

       if (length($row[$std_columns[$k]]) > 0)
       {
	  $not_all_empty = 1;
	  push(@tmp_list,$row[$std_columns[$k]]);
       }
    }

    if ($not_all_empty == 1)
    {
       my $first = 0;
       my $last = @tmp_list - 1;
       $first = ($first > $last) ? $last : $first;

       if (($last >= $first) and ($last < @tmp_list))
       {
	  $str = 0;
	  my $sum_x = 0;
	  my $sum_diffs = "";
	  my $avg_s = 0;
	  my $entries = $last - $first +1;
	  for (my $k = $first; $k <= $last; $k++)
	  {
	     $sum_x = $sum_x + $tmp_list[$k];
	  }
	  my $avg_s = ($entries > 0) ? $sum_x / $entries : "";;
	  if (length($avg_s) > 0)
	  {
	      for (my $k = $first; $k <= $last; $k++)
	      {
		  $sum_diffs = $sum_diffs + ($tmp_list[$k] - $avg_s)*($tmp_list[$k] - $avg_s);
	      }
	  }
	  $str = (($entries > 1) and (length($sum_diffs) > 0)) ? sqrt($sum_diffs / ($entries-1)) : 0;
       }
    }
    $str = (($not_all_empty == 1) and (length($str) > 0)) ? &format_number($str, $significant_nums) : "";
  }
elsif (length($stderr_columns_string) > 0)
  {
    my @row = split(/\t/);

    $str = "";
    my $not_all_empty = 0;
    my @tmp_list = ();

    for(my $k=0; $k < @stderr_columns; $k++)
    {
       next if $ignore_zero_entries and $row[$stderr_columns[$k]] == 0;

       if (length($row[$stderr_columns[$k]]) > 0)
       {
	  $not_all_empty = 1;
	  push(@tmp_list,$row[$stderr_columns[$k]]);
       }
    }

    if ($not_all_empty == 1)
    {
       my $first = 0;
       my $last = @tmp_list - 1;
       $first = ($first > $last) ? $last : $first;

       if (($last >= $first) and ($last < @tmp_list))
       {
	  $str = 0;
	  my $sum_x = 0;
	  my $sum_diffs = "";
	  my $avg_s = 0;
	  my $entries = $last - $first +1;
	  for (my $k = $first; $k <= $last; $k++)
	  {
	     $sum_x = $sum_x + $tmp_list[$k];
	  }
	  my $avg_s = ($entries > 0) ? $sum_x / $entries : "";;
	  if (length($avg_s) > 0)
	  {
	      for (my $k = $first; $k <= $last; $k++)
	      {
		  $sum_diffs = $sum_diffs + ($tmp_list[$k] - $avg_s)*($tmp_list[$k] - $avg_s);
	      }
	  }
	  $str = (($entries > 1) and (length($sum_diffs) > 0)) ? sqrt($sum_diffs / ($entries*($entries-1))) : 0;
       }
    }
    $str = (($not_all_empty == 1) and (length($str) > 0)) ? &format_number($str, $significant_nums) : "";
  }
  elsif ($start_counter == 1) { $str = $row_counter; }

  if ($column_location != -1)
  {
    my @row = split(/\t/);

    for (my $i = 0; $i < @row; $i++)
    {
      if ($i == $column_location)
      {
	if ($i == 0) { print "$str\t"; }
	else { print "\t$str"; }
      }

      if ($i > 0) { print "\t"; }

      print "$row[$i]";
    }
    print "\n";
  }
  elsif ($beginning == 1) { print "$str\t$_\n"; }
  else { print "$_\t$str\n"; }

  $row_counter++;
}

__DATA__

add_column.pl <file>

   Adds a column to each of the lines of a file.
   Operations on more than two columns can be specified using the ',-' notation (e.g., 1,4-6,23 means columns 1,4,5,6,23)

   -b:                 add the column as the first in the file (default: add at the end)
   -c <num>:           add the column before column number <num>
   -e:                 add an empty column.
   -n <num>            add a column counter, starting at num.
   -s <str>            add a column with the specified string
   -d <c1,c2>          add a column which is the value of column1 / column2
   -u <c1,c2>          add a column which is the value of column1 - column2
   -cat <c1,c2,..,ck>  add a column which is the value of column1 . column2 . ... . columnk
   -delim <str>        to use with cat option. Concatenate columns with delim between them
   -m <c1,c2,..,ck>    add a column which is the value of column1 * column2 *..* columnk
   -a <c1,c2,..,ck>    add a column which is the value of column1 + column2 +..+ columnk
   -min <c1,c2,..,ck>  add a column which is the value of min(column1,column2,..,columnk)
   -max <c1,c2,..,ck>  add a column which is the value of max(column1,column2,..,columnk)
   -count <c1,c2,..,ck>add a column which is the count of non empty entries over column1,column2,..,columnk.
   -ave <c1,c2,..,ck>  add a column which is the value of average(column1,column2,..,columnk) !!
   -quant <low,high>   the quantile range of values on which to perform the operation (default: 0,1)
                       !! currently works only with -ave !!
                       e.g., to compute a trimmed mean in quantiles range 0.1 to 0.8 use: -ave c1,..,ck -quant 0.1,0.8
   -std <c1,c2,..,ck>  add a column which is the value of std(column1,column2,..,columnk) 
   -se <c1,c2,..,ck>   add a column which is the value of std(column1,column2,..,columnk) 

   -nozero:            use only non-zero column values (out of the given columns list).

   -f <name>:          add the column from the specified file
   -sn <int>:          The significant numbers to print (default: 2).

