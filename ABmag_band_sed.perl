#!/usr/bin/perl -w
# simple script to call src_ray with approximate arguments to yield the right number of input photons for psf shape etc..

use POSIX;

($band,$alpha,$mag,$exptime)=("i",2,26,15);
$extra_args="";

$slop=0.03;

while (($#ARGV >0) && ($_=$ARGV[0],$_=~/^-/)) {
    shift;
    /-b/ && ($band=$ARGV[0],shift,next);
    /-p/ && ($alpha=$ARGV[0],shift,next);
    /-m/ && ($mag=$ARGV[0],shift,next);
    /-t/ && ($exptime=$ARGV[0],shift,next);
    /-e/ && ($extra_args.=" -e",next);
    /-s/ && ($slop=$ARGV[0],shift,next);
    /-S/ && ($seeing=$ARGV[0],shift,$extra_args.=" -S $seeing",next);
    printf STDERR "unknown switch: $ARGV[0]\n";
    exit(1);
}


if ($band =~ /^[uU]/) {
    $band="U";
    $slop=0.0;
    @blims=(3210/(1+$slop),3915*(1+$slop));
    $bcen=sqrt($blims[$[]*$blims[$[+1]);
    $slop=0.03;
    @blims=(3210/(1+$slop),3915*(1+$slop));
} elsif ($band =~ /^[gG]/) {
    $band="G"; 
    $slop=0.0;
    @blims=(3997/(1+$slop),5500*(1+$slop));
    $bcen=sqrt($blims[$[]*$blims[$[+1]);
    $slop=0.03;
    @blims=(3997/(1+$slop),5500*(1+$slop));
} elsif ($band =~ /^[rR]/) {
    $band="R"; 
    $slop=0.0;
    @blims=(5490/(1+$slop),6900*(1+$slop));
    $bcen=sqrt($blims[$[]*$blims[$[+1]);
    $slop=0.03;
    @blims=(5490/(1+$slop),6900*(1+$slop));
} elsif ($band =~ /^[iI]/) {
    $band="I"; 
    $slop=0.0;
    @blims=(6880/(1+$slop),8200*(1+$slop));
    $bcen=sqrt($blims[$[]*$blims[$[+1]);
    $slop=0.03;
    @blims=(6880/(1+$slop),8200*(1+$slop));
} elsif ($band =~ /^[zZ]/) {
    $band="Z"; 
    $slop=0.0;
    @blims=(8185/(1+$slop),9225*(1+$slop));
    $bcen=sqrt($blims[$[]*$blims[$[+1]);
    $slop=0.03;
    @blims=(8185/(1+$slop),9225*(1+$slop));
} elsif ($band =~ /^[yY]2/) {
    $band="Y2"; 
    $slop=0.0;
    @blims=(9790/(1+$slop),10245*(1+$slop));
    $bcen=sqrt($blims[$[]*$blims[$[+1]);
    $slop=0.03;
    @blims=(9790/(1+$slop),10245*(1+$slop));
} elsif ($band =~ /^[yY]3/) {
    $band="Y3"; 
    $slop=0.0;
    @blims=(9745/(1+$slop),11220*(1+$slop));
    $bcen=sqrt($blims[$[]*$blims[$[+1]);
    $slop=0.03;
    @blims=(9745/(1+$slop),11220*(1+$slop));
} elsif ($band =~ /^[yY]4/) {
    $band="Y4"; 
    $slop=0.0;
    @blims=(9225/(1+$slop),12000*(1+$slop));
    $bcen=sqrt($blims[$[]*$blims[$[+1]);
    $slop=0.03;
    @blims=(9225/(1+$slop),12000*(1+$slop));
} else {
    printf STDERR "bad band choice: $band choose from [ugrizy{2,3,4}]\n";
    exit(1);
}

# exp(log(10)*arg) = 10^(arg)
$norm=exp(log(10)*(-7.36022-0.4*($mag-26)-log($bcen/5000)/log(10)));
printf STDERR "norm = $norm\n";

# exit;

$norm=$norm*(pow(12.398,$alpha-1)*pow($bcen,2-$alpha));

$cmd="src_ray $extra_args -b $blims[$[] $blims[$[+1] $alpha -f $norm -r 4500 -t $exptime -s 0 0";
$cmd .= " -S $seeing" if (defined($seeing) && $seeing>0);
printf STDERR "band $band alpha $alpha norm $norm\n";
printf STDERR "$cmd\n";

open(F,"$cmd |") || die;;
print while (<F>);
close(F);
# done.
