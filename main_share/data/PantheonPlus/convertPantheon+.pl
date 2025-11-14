#!/usr/bin/perl

my $data_table = "Pantheon+SH0ES.dat";
my $cov_matrix = "Pantheon+SH0ES_STAT+SYS.cov";

my @sigma;

sub convertTable
{
  open (TABLE, "<$data_table") || die "cannot open $data_table: $!\n";
  open (OUTPUT, ">pantheon+Sample.d");
  while (<TABLE>)
  {
    if (!/^#/)
    {
      my @line = split;
      push (@sigma, $line[11]);
      printf OUTPUT
        ("%6.4f  %8.4f  %8.4f  %8.4f\n",
         $line[4], $line[10], $line[11]/2, $line[11]/2);
    }
  }
  close (OUTPUT);
  close (TABLE);
}

sub convertMatrix
{
  open (MATRIX, "<$cov_matrix") || die "cannot open $cov_matrix: $!\n";
  open (OUTPUT, ">pantheon+Sample_CM.d");
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
