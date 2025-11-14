#!/usr/bin/perl

# Perl script for converting the csv table of BAO measurements to a TeX
# table or to a data table

sub printTeXHeader {
  print TEXF <<EOH;
\\begin{table*}
  \\caption{BAO data}
  \\label{tab:2a}
  \\begin{tabular}{llllp{70mm}l}
    \$n\$ & \$z\$ & \$D_\\mathrm{A}/r_\\mathrm{d}\$ & \$\\Delta(D_\\mathrm{A}/r_\\mathrm{d})\$ &
    Description & Reference \\\\
    \\hline
EOH
}

sub printTeXFooter {
print TEXF <<EOF
  \\end{tabular}
\\end{table*}
EOF
}

sub printTblHeader {
print TBLF <<EOH;
# Cleaned data set from Sofía Vidal's Bachelor thesis
#
# r_d: sound horizon at decoupling
# column 1: redshift
# column 2: distance modulus from angular diameter distance (D_A/r_d)
# column 3: uncertainty in distance modulus
#
EOH
}

sub csv2TeX {
  my $nline = shift;
  my $tmpfile = "/tmp/_sorted.csv";
  system ("sort -n -t';' +2 -o $tmpfile baoTable.csv");

  my $format='    %5.3f & %9.6f & %8.6f &\n    %s &\n    \\cite{%s} \\\\\n';
  my $count = 0;
  my $line = 0;

  open (FILE, "<$tmpfile");
  open (TEXF, ">baoTable.tex");

  printTeXHeader;
  while (<FILE>) {
    if (!/ignore/) {
      chomp;
      my @line = split (';');
      $count++;
      $line++;
      $line[6] =~ s/A&A/A+A/;
      printf TEXF
        ("    %2i & %5.3f & %7.4f & %6.4f &\n    %s &\n    \\cite{%s} \\\\\n",
         $count, $line[2], $line[3], $line[4], $line[5], $line[6]);
      if ($line == $nline) {
        $line = 0;
        printTeXFooter;
        printTeXHeader;
      }
    }
  }
  printTeXFooter;

  close (TEXF);
  close (FILE);

  system ("rm -f $tmpfile");
}

sub csv2tbl {
  open (FILE, "<baoTable.csv");
  open (TBLF, ">baoSample.d");
  
  printTblHeader;
  my $factor = 5.0/log (10.0);
  while (<FILE>) {
    if (!/ignore/) {
      chomp;
      my @line = split (';');
      # convert distance to distance modulus
      my $mu = $factor*log ((1.0+$line[2])**2*$line[3]*1.0e5);
      # convert distance uncertainties to uncertainties in distance modulus
      my $delta_mu = 0.5*$factor*$line[4]/$line[3];
      printf TBLF
        ("%6.4f  %7.4f  %6.4f  %6.4f\n", $line[2], $mu, $delta_mu, $delta_mu);
    }
  }
  
  close (TBLF);
  close (FILE);
}

sub addOffset {
  my $findOffset = "./findOffset";
  my $clean = 0;
  if (! -x $findOffset.".x") {
    if (! -r $findOffset.".cpp") {
      die "offset code not found, exiting\n";
    } else {
      my $command = `grep g++ $findOffset".cpp" | sed s/^.*g++/g++/;`;
      system ($command);
      $clean = 1;
    }
  }
  my $flag = 1;
  my $tmpfile = "/tmp/_baoSample.d";

  open (FILE, "<baoSample.d");
  open (OUTP, ">$tmpfile");
  while (<FILE>) {
    if (/^#/) { print OUTP; }
    else {
      if ($flag) {
        $flag = 0;
        my $entry = `$findOffset".x"`;
        print OUTP $entry;
      }
      print OUTP;
    }
  }
  close (FILE);
  close (OUTP);
  system ("mv $tmpfile baoSample.d");
  if ($clean) { system ("rm $findOffset.x"); }
}

# ----------

my $nline = 25;
my $texflag = 0;

while (@ARGV) {
  my $arg = shift (@ARGV);
  if ($arg eq "-n") { $nline = shift (@ARGV); }
  if ($arg == "-t") { $texflag = 1; }
}

if ($texflag) {
  csv2TeX ($nline);
} else {
  csv2tbl;
  addOffset;
}

exit 0;
