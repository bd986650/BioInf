version 1.0

workflow Hello {
  call SayHello
}

task SayHello {
  command <<<
    echo "Hello WDL"
  >>>

  output {
    String msg = read_string(stdout())
  }
}