#!/usr/bin/perl
#---------------------------------------- V3.1 (C) Ench, Nov 2005, Jan 2013 --
#
## sumablastplus - blast+ output summarization
#
#| Usage: sumablast [-opts] inpfn.txt outfn.tsv
#| Options:
#|   -maxhits n   - maximum number of hits per query to include in output
#|   -empties     - include empty queries (with no hits found)
#|   inpfn.txt    - input file (text, should be output from blast program)
#|   outfn.tsv    - output file (tab separated)
#|
#-----------------------------------------------------------------------------
use strict;

my $DEBUG=($ENV{DEBUG}||0);
my $QUERYRE="Query= ?";
my @CSVHDR=qw(Query_ID Query_length Target_ID Target_length
  Target_description E-value Identity% Aligned_seq_length Score
  Positives% Gaps%);
my($val,$emptyq,$maxhits,$inpfn,$outfn);
my($INP,$OUT,$nextln,@ln,$maxcol,$hit);
my($qid,$qlen,$tid,$tlen,$tdesc,$eval,$idp,$alen,$score,$positives,$gaps);

usage() if (!@ARGV);
while (defined($_=shift)) {s/^--/-/; 
  if (/^-.+=/) {($_,$val)=split('=');$val=unshift(@ARGV,$val)}
  else {$val=(defined($ARGV[0]) and $ARGV[0]!~/^-/)?1:undef}
  if (opt($_,"-help",2)) {usage()}
  elsif (opt($_,"-empties",2)) {$emptyq=1}
  elsif (opt($_,"-maxhits",2)) {$maxhits=($val?shift:" ")}
  elsif (/^-.+/) {err("Illegal option ($_)")}
  elsif (!defined($inpfn)) {$inpfn=$_}
  elsif (!defined($outfn)) {$outfn=$_}
  else {usage()}
}
$maxhits=5 if (!defined($maxhits));
err("illegal value (-maxhits $maxhits)") if ($maxhits=~/\D/ or $maxhits<1);

if (!defined($inpfn) or !defined($outfn)) {usage()}
elsif ($inpfn eq "-") {$INP=\*STDIN}
elsif (!open($INP,"<",$inpfn)) {err("$inpfn: $!")}
if ($outfn eq "-") {$OUT=\*STDOUT}
elsif ($outfn!~/^[\w\.\-\/]+$/) {err("Illegal filename ($outfn)")}
elsif (!open($OUT,">",$outfn)) {err("$outfn: $!")}

readheader();
while (readquery()) {
# debug("READQUERY:'%s'",$ln[0]);
  $qid=scan4re($QUERYRE,"Length=");
  ($qlen=shift(@ln))=~s/Length= *(\d+)/$1/;
# debug("QUERY_ID=%s",$qid); debug("QUERY_LEN=%d",$qlen);
  if (!skip2re(">")) {wrrow($qid,$qlen) if ($emptyq); next}
  $hit=0;
  while (@ln) {$hit++;
    if ($hit>$maxhits) {skip2re(">");next} # debug(" HIT=%d (SKIPPED)",$hit);
#   debug(" HIT=%d",$hit);
    $_=scan4re(">","Length="); s/\s+/ /g;
    ($tid,$tdesc)=split(' ',$_,2);
    ($tlen=shift(@ln))=~s/Length= *(\d+)/$1/;
#   debug(" TARGET_ID=%s",$tid);
#   debug(" TARGET_DESC=%s",$tdesc);
#   debug(" TARGET_LEN=%s",$tlen);
    $_=scan4re(" ","Query "); s/\s+/ /g;
#   debug(" SCORELINE='%s'",$_);
    ($score)=(/Score = +([\d\.]+) bits /);
    ($eval)=(/Expect(?:\(\d+\))? = ([^\s,]+)/);
    ($alen,$idp)=(/ Identities = \d+\/(\d+) \((\d+)%\)/);
    ($positives)=(/ Positives = \d+\/\d+ \((\d+)%\)/);
    ($gaps)=(/ Gaps = \d+\/\d+ \((\d+)%\)/);
#   debug("  EXPVAL=%s",$eval);
#   debug("  IDENT%%=%s",($idp||""));
#   debug("  ALIGNLEN=%s",($alen||""));
#   debug("  SCORE=%s",$score);
#   debug("  POSITIVES%%=%s",$positives);
#   debug("  GAPS%%=%s",$gaps);
    wrrow($qid,$qlen,$tid,$tlen,$tdesc,$eval,$idp,$alen,$score,$positives,$gaps);
    skip2re(">");
  }
}
close($OUT);
close($INP);
#-----------------------------------------------------------------------------
sub wrrow { my(@row)=@_; my($i);
  if (!$maxcol) {$maxcol=@CSVHDR; print $OUT join("\t","",@CSVHDR),"\n"}
  for ($i=0;$i<$maxcol;$i++) {$row[$i]="" if (!defined($row[$i]))}
  print $OUT join("\t","",@row),"\n";
# debug("WRROW:<%s>",join("><",@row));
}
#-----------------------------------------------------------------------------
sub scan4re { my($startre,$endre)=@_; my($s);
# debug2("SCAN4RE(/%s/,/%s/) >(%s)",$startre,$endre,$ln[0]);
  err("Unable to find '$startre' ($ln[0])") if ($ln[0]!~/^$startre/);
  ($s=shift(@ln))=~s/^$startre *//;
  chomp($s.=" ".shift(@ln)) while (@ln and $ln[0]!~/^$endre/);
  err("End of data while looking for /$endre/") if (!@ln);
# debug2("-->(%s)",$s);
  $s;
}
#-----------------------------------------------------------------------------
sub skip2re { my($re,$fatal)=@_;
# debug2("SKIP2RE(/%s/) >(%s)",$re,$ln[0]);
  shift(@ln);
  while (@ln and $ln[0]!~/^$re/) {shift(@ln)} # debug2(" (%s)",$ln[0]);
  err("End of data while looking for /$re/") if ($fatal and !@ln);
# debug2("-->(%s)",($ln[0]||""));
  scalar(@ln);
}
#-----------------------------------------------------------------------------
sub readheader { my($magic)='T?BLAST[PNX] \d+\.\d+\.\d+\+';
  err("$inpfn: Empty file") if (!defined($_=<$INP>));
  chomp(),err("Wrong magic (".substr($_,0,30).")") if ($_!~/^$magic/);
  while (defined($_=<$INP>) and !/^$QUERYRE/) {}
  err("Unable to find query") if (!$_);
  chomp($nextln=$_);
}
#-----------------------------------------------------------------------------
sub readquery { @ln=();
  push @ln,$nextln if ($nextln); $nextln=undef;
  while (<$INP>) {chomp; s/\r$//; s/\s/ /g; next if (/^\s*$/);
    $nextln=$_,last if (/^$QUERYRE/);
    push @ln,$_;
  }
  scalar(@ln);
}
#-----------------------------------------------------------------------------
sub debug { my($f,@p)=@_; printf("|DBG|$f\n",@p) if ($DEBUG)}
sub debug2 { my($f,@p)=@_; printf("|DBG|$f\n",@p) if ($DEBUG and $DEBUG>1)}
sub err {my($s)=@_; die sprintf("%s: %s\n",(split('/',$0))[-1],$s)}
sub opt {my($s,$k,$m)=@_; length($s)>=($m?$m:length($k)) and $k=~/^\Q$s\E/}
sub usage {open(FN,"<$0") && print grep {s/^#\| ?//} <FN>; close(FN); exit(1)}
sub unquote {my($s)=@_; $s=~s/^(['"])(.*)\1$/$2/; return($s)}
sub dirname {my($p)=@_; $p=~s|[^/]+/?$||; $p=~s|/$|| if ($p ne "/");$p?$p:"."}
#-----------------------------------------------------------------------------
#sub DUMPHASH {my %h=@_;printf "{%s}='%s'\n",$_,$h{$_} foreach (sort keys %h)}
#sub DUMPARR {for (my $i=0;$i<@_;$i++) {printf "[%d]='%s'\n",$i,$_[$i]}}
#sub ECHO { my %c=qw(W 30;47 R 37;41 G 30;42 Y 30;43 B 37;44 M 37;45 C 30;46);
#  my($s)=@_; my($n)=($s=~s/\+$//); $s=~s/^([a-z]?)(\d?)://i;
#  printf "%*s\e[%sm%s\e[m".($n?"":"\n"),$2?$2:0,"",$c{uc($1?$1:"C")},$s;
#}
#-----------------------------------------------------------------------------
