#include <stdlib.h>

#include <End.hh>
#include <RAT/AnyParse.hh>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
  auto parser = new RAT::AnyParse(argc, argv);
  std::cout << "End version: " << RAT::ENDVERSION << std::endl;
  auto end = END::End(parser, argc, argv);
  end.Begin();
  end.Report();
}
