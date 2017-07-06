#include "CLI/CLI.hpp"
#include "FastxParser.hpp"
#include "string_view.hpp"
#include "FileOps.hpp"
#include "ScopedTimer.hpp"
#include <cmath>
#include <iterator>
#include <iostream>
#include <vector>


void createFastqFiles(std::string fastaFile,
                      std::string outdir,
                      bool isPairedEnd,
                      uint32_t readLen){
  {
    ScopedTimer st ;
    std::vector<std::string> read_file = {fastaFile} ;
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(read_file, 1, 1) ;
    parser.start() ;


    if(outdir.back() == '/')
      outdir.pop_back() ;

    MakeDir(outdir.c_str()) ;
    
    std::ofstream fastqFile, fastqFile_1, fastqFile_2;
    if(!isPairedEnd){
      fastqFile.open(outdir + "/1.fq") ;
    }else{
      fastqFile_1.open(outdir + "/1.fq") ;
      fastqFile_2.open(outdir + "/2.fq") ;
    }

   

    size_t globalId{0};
    size_t rn{0} ;
    auto rg = parser.getReadGroup() ;
    while(parser.refill(rg)){
      for(auto& rp: rg){

        if(rn % 5000 == 0){
          std::cerr << "rn: "<<rn << "\n" ;
        }
        ++rn ;
        auto& r1 = rp.seq ;
        auto& name = rp.name ;

        stx::string_view r1view(r1) ;

        if(r1.length() > readLen){
          if(!isPairedEnd){
            size_t index{0} ;
            for(;index < (r1view.size()-readLen+1);++index){
              fastqFile << "@"<<name<<"_"<<globalId<<"\n" ;
              fastqFile << r1view.substr(index,readLen) << "\n";
              fastqFile << "+" << "\n" ;
              std::string quality(readLen,'I') ;
              fastqFile << quality << "\n" ;
              ++globalId ;
            }
          }else{
            size_t index{0} ;
            for(;index < (r1view.size()-readLen);++index){
              fastqFile_1 <<"@"<<name<<"_"<<globalId<<"/1\n" ;
              fastqFile_2 <<"@"<<name<<"_"<<globalId<<"/2\n" ;
              fastqFile_1 << r1view.substr(index,readLen) << "\n";
              fastqFile_2 << r1view.substr(index+1,readLen) << "\n";
              fastqFile_1 << "+" << "\n" ;
              fastqFile_2 << "+" << "\n" ;
              std::string quality(readLen,'I') ;
              fastqFile_1 << quality << "\n" ;
              fastqFile_2 << quality << "\n" ;
              ++globalId ;
            }
            
          }

        }else{

          std::string padstr(readLen - r1.length(),'N') ;
          std::string readStr = r1+padstr ;

          if(!isPairedEnd){
            fastqFile <<"@"<<name<<"_"<<globalId<<"\n" ;
              fastqFile << readStr << "\n";
              fastqFile << "+" << "\n" ;
              std::string quality(readLen,'I') ;
              fastqFile << quality << "\n" ;
              ++globalId ;
          }else{
            fastqFile_1 <<"@"<<name<<"_"<<globalId<<"/1\n" ;
            fastqFile_2 <<"@"<<name<<"_"<<globalId<<"/2\n" ;
              fastqFile_1 << readStr << "\n";
              fastqFile_2 << readStr << "\n";
              fastqFile_1 << "+" << "\n" ;
              fastqFile_2 << "+" << "\n" ;
              std::string quality(readLen,'I') ;
              fastqFile_1 << quality << "\n" ;
              fastqFile_2 << quality << "\n" ;
              ++globalId ;
          }
        }

      }
    }

    if(!isPairedEnd){
      fastqFile.close() ;
    }else{
      fastqFile_1.close() ;
      fastqFile_2.close() ;
    }
  }
}


int main(int argc, char* argv[]){
  CLI::App app{"Playing tools"} ;
  auto fastaqApp = app.add_subcommand("fastqgen", "fastq generator") ;

  //options for creating the required files
  uint32_t readLen{31} ;
  std::string outdir ;
  std::string fastaFile ;
  bool isPairedEnd{false} ;
  //end

  fastaqApp
    ->add_option("-f,--fasta", fastaFile,"the full path to the fasta file")
    ->required() ;

  fastaqApp
    ->add_option("-o,--outdir", outdir,"output directory")
    ->required() ;

  fastaqApp
    ->add_flag("-p,--paireend", isPairedEnd,
               "library type");
  //fastaqApp
  //->add_flag("-p,--library", isPairedEnd, "paired end or not") ;

  fastaqApp
    ->add_option("-r,--read-len", readLen,"read length") ;

  try{
    app.parse(argc,argv) ;
  } catch (const CLI::ParseError& e){
    std::cerr << "\n Not proper arguments\n " ;
    app.exit(e) ;
  }

  if(app.got_subcommand(fastaqApp)){
    createFastqFiles(fastaFile, outdir, isPairedEnd, readLen) ;
  }else{
    std::cerr << "wrong subcommand\n" ;
  }


}
