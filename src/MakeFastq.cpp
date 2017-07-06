#include "CLI/CLI.hpp"
#include "FastxParser.hpp"
#include "string_view.hpp"
#include "FileOps.hpp"
#include "ScopedTimer.hpp"
#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
#include <cmath>
#include <iterator>
#include <iostream>
#include <vector>


void createFastqFiles(std::string fastaFile,
                      std::string outdir,
                      bool isPairedEnd,
                      uint32_t readLen,
                      bool ignoreShortSeq){
  
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

    CanonicalKmer::k(readLen) ;
    pufferfish::CanonicalKmerIterator kit_end ;

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
          pufferfish::CanonicalKmerIterator kit1(r1) ;
          if(!isPairedEnd){

            for(; kit1 != kit_end ; ++kit1){
              fastqFile << "@"<<name<<"_"<<globalId<<"\n" ;
              fastqFile << kit1->first.to_str() << "\n";
              fastqFile << "+" << "\n" ;
              std::string quality(readLen,'I') ;
              fastqFile << quality << "\n" ;
              ++globalId ;
            }

          }else{
            for(; kit1 != kit_end ; ){
              auto kit2 = kit1 ;
              fastqFile << "@"<<name<<"_"<<globalId<<"\n" ;
              fastqFile_1 <<"@"<<name<<"_"<<globalId<<"/1\n" ;
              fastqFile_2 <<"@"<<name<<"_"<<globalId<<"/2\n" ;
              fastqFile_1 << kit1->first.to_str() << "\n";
              ++kit1;
              if(kit1 != kit_end){
                fastqFile_2 << kit1->first.to_str() << "\n";
                ++kit1 ;
              }else{
                fastqFile_2 << kit2->first.to_str() << "\n";
                break ;
              }
              fastqFile_1 << "+" << "\n" ;
              fastqFile_2 << "+" << "\n" ;
              std::string quality(readLen,'I') ;
              fastqFile_1 << quality << "\n" ;
              fastqFile_2 << quality << "\n" ;
              ++globalId ;
            }
            
          }

        }else{
          if(ignoreShortSeq)
            continue ;

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
  bool ignoreShortSeq{true} ;
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

  fastaqApp
    ->add_flag("-i,--ignore", ignoreShortSeq,
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
    createFastqFiles(fastaFile, outdir, isPairedEnd, readLen, ignoreShortSeq) ;
  }else{
    std::cerr << "wrong subcommand\n" ;
  }


}
