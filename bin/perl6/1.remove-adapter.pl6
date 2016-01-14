sub usage(Bool :$verbose=False){
  say "Hello, this is the usage page";
}
sub MAIN(Str :$adapter!, int :$min-length=12, int :$max-length=30, int :$min-quality=20, Str :$prefix!, Str :$outdir!, Bool :$help){
  if $help {
    &usage();
    fail "examples";
  }
  if $adapter > 15{
    my $adapter-new = $adapter.substr(0,15);
  }
}
