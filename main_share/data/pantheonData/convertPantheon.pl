#!/usr/bin/perl

my $data_table = "hlsp_ps1cosmo_panstarrs_gpc1_all_model_v1_lcparam-full.txt";
my $cov_matrix = "hlsp_ps1cosmo_panstarrs_gpc1_all_model_v1_sys-full.txt";

my @sigma;

sub convertTable
{
  my $root2 = sqrt (2.0);
  open (TABLE, "<$data_table") || die "cannot open $data_table: $!\n";
  open (OUTPUT, ">pantheonSample.d");
  while (<TABLE>)
  {
    if (!/^#/)
    {
      my @line = split;
      push (@sigma, $line[5]);
      printf OUTPUT
        ("%6.4f  %8.4f  %8.4f  %8.4f\n",
         $line[1], $line[4], $line[5]/2, $line[5]/2);
    }
  }
  close (OUTPUT);
  close (TABLE);
}

sub convertMatrix
{
  open (MATRIX, "<$cov_matrix") || die "cannot open $cov_matrix: $!\n";
  open (OUTPUT, ">pantheonSample_CM.d");
  my $i = 0;
  my $j = 0;
  my $n = 0;
  while (<MATRIX>)
  {
    chomp;
    my @line = split;
    if ($i == 0)
      { $n = $line[0]; }
    else
    {
      my $cov = $line[0];
      if (($i-1) % $n == $j)
        { $cov += $sigma[$j]*$sigma[$j]; }
      print OUTPUT "$cov  ";
    }
    $i++;
    if ($i > 1 && ($i-1) % $n == 0)
    {
      $j++;
      print OUTPUT "\n";
    }
  }
  close (OUTPUT);
  close (MATRIX);
}

convertTable;
convertMatrix;

exit (0);
