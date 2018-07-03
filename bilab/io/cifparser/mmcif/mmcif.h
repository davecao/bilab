//
//  cif.h
//  biTKCXX
//
//  Created by 曹巍 on 2018/06/19.
//  Copyright © 2018年 巍 曹. All rights reserved.
//
// This idea was inspired by cif.hpp developped Global Phasing Ltd
//
// Copyright 2017 Global Phasing Ltd.
//
// CIF parser (based on PEGTL), struct Document that represents the CIF file,
// and a set of actions for the parser to prepare Document.

#ifndef mmcif_h
#define mmcif_h


#include "tao/pegtl.hpp"
//#include <tao/pegtl/contrib/tracer.hpp>  // for debugging


#include "mmcif_rules.h"
#include "CIFDocument.h"
#include "utils.h"

namespace bilab{
  namespace pegtl = tao::pegtl;
  namespace mmcif{
    using namespace pegtl;
    //namespace rules = bilab::mmcif::rules;
    // error message
    template<typename Rule>
    struct Errors
    : public pegtl::normal<Rule>
    {
      static const std::string msg;
      template<typename Input, typename ... States>
      static void raise(const Input& in, States&& ...) {
        throw pegtl::parse_error(msg, in);
      }
    };
    
#if defined(_MSC_VER)
# define error_msg(x) \
template<> __declspec(selectany) const std::string Errors<x>::msg
#else
# define error_msg(x) \
template<> const std::string Errors<x>::msg __attribute__((weak))
#endif
  error_msg(rules::quoted_tail<one<'\''>>) = "unterminated 'string'";
  error_msg(rules::quoted_tail<one<'"'>>) = "unterminated \"string\"";
  error_msg(pegtl::until<rules::field_sep>) = "unterminated text field";
  error_msg(rules::value) = "expected value";
  error_msg(rules::datablockname) = "unnamed DataBlock";
  error_msg(rules::framename) = "unnamed save_ frame";
#undef error_msg
    template<typename T>
    const std::string Errors<T>::msg = "parse error";
    
    /**** parsing actions that fill the storage ****/
    template<typename Rule>
    struct Action : pegtl::nothing<Rule>
    {};
 
    template<>
    struct Action<rules::datablockname> {
      template<typename Input>
      static void apply(const Input& in, CIFDocument& out) {
        //std::cout<< "datablockname: " << in.string() <<std::endl;
        out.source = in.string();
        //out.blocks.emplace_back(in.string());
        //out.items_ = &out.blocks.back().items;
      }
    };
    
    template<>
    struct Action<rules::framename> {
      template<typename Input>
      static void apply(const Input& in, CIFDocument& out) {
        //out.line_number += in.iterator().line;
        //std::cout << "framename: " << in.string() <<"\n";
        //out.items_->emplace_back(FrameArg{in.string()});
        //out.items_->back().line_number = in.iterator().line;
        //out.items_ = &out.items_->back().frame.items;
      }
    };
    
    template<>
    struct Action<rules::endframe> {
      template<typename Input>
      static void apply(const Input& in, CIFDocument& out) {
        //out.line_number += in.iterator().line;
        //std::cout << "endframe: " << in.string() <<"\n";
        //out.items_ = &out.blocks.back().items;
      }
    };
    
    template<>
    struct Action<rules::str_global> {
      template<typename Input>
      static void apply(const Input& in, CIFDocument& out) {
        
        //out.line_number += in.iterator().line;
        //std::cout << "str_global: " << in.string() <<"\n";
        //out.blocks.emplace_back();
        //out.items_ = &out.blocks.back().items;
      }
    };
    
    template<>
    struct Action<rules::tag> {
      template<typename Input>
      static void apply(const Input& in, CIFDocument& out) {
        out.currTag = in.string();
        out.line_number += in.iterator().line;
        out._tags.emplace_back(in.string());
      }
    };
    
    template<>
    struct Action<rules::value> {
      template<typename Input>
      static void apply(const Input& in, CIFDocument& out) {
        //std::cout << "(tag, value): " << out.currTag <<","<< in.string() <<"\n";
        out._values.emplace_back(in.string());
      }
    };
    
    template<>
    struct Action<rules::str_loop> {
      template<typename Input>
      static void apply(const Input& in, CIFDocument& out) {
        // Clear _tags and _values
        out._values.clear();
        out._tags.clear();
        out.isloop = true;
        out.line_number = static_cast<int>(in.iterator().line);
      }
    };
    
    template<>
    struct Action<rules::loop_tag> {
      template<typename Input>
      static void apply(const Input& in, CIFDocument& out) {
        out._tags.emplace_back(in.string());
        out.currTag = in.string();
      }
    };
    
    template<>
    struct Action<rules::loop_value> {
      template<typename Input>
      static void apply(const Input& in,CIFDocument& out) {
        out._values.emplace_back(in.string());
      }
    };
    
    template<>
    struct Action<rules::loop> {
      template<typename Input>
      static void apply(const Input& in, CIFDocument& out) {
        // Here, all the information in the block
        // e.g.,
        //    loop_
        //    _database_2.database_id
        //    _database_2.database_code
        //    PDB   1A0S
        //    WWPDB D_1000170250
        //#
      }
    };
    
    template<>
    struct Action<rules::comment> {
      template<typename Input>
      static void apply(const Input& in, CIFDocument& out) {
        // comment
        //if (!(out._values.empty() && out._tags.empty()) && out.isloop) {
          out.add_block_loop();
          out.isloop = false;
          //}
      }
    };
    
    template<typename Input>
    bool parse_input(CIFDocument& d, Input&& in) {
      // Call pegtl's parse
      bool state = pegtl::parse<rules::file, Action, Errors>(in, d);
      d.source = in.source();
      //check_duplicates(d);
      return state;
    }
    
    template<typename Input>
    bool read_input(CIFDocument& d, Input&& in) {
      //CIFDocument* doc = new CIFDocument();
      bool state = parse_input(d, in);
      //return doc;
      return state;
    }
    
    inline bool read_file(const std::string& filename, CIFDocument& d) {
      pegtl::file_input<> in(filename);
      bool state = read_input(d, in);
      return state;
    }
    
    inline bool read_string(const std::string& data, CIFDocument& d) {
      pegtl::memory_input<> in(data, "string");
      bool state = read_input(d, in);
      return state;
    }
    
    inline bool read_memory(const char* data, size_t size, const char* name, CIFDocument& d) {
      pegtl::memory_input<> in(data, size, name);
      bool state = read_input(d, in);
      return state;
    }
    
   inline bool read_cstream(std::FILE *f, size_t maximum, const char* name, CIFDocument& d) {
      pegtl::cstream_input<> in(f, maximum, name);
      bool state = read_input(d, in);
     return state;
    }
    
    inline bool read_istream(std::istream &is,
                      size_t maximum,
                      const char* name,
                      CIFDocument& d) {
      pegtl::istream_input<> in(is, maximum, name);
      bool state = read_input(d, in);
      return state;
    }
    /*
    template<typename Input>
    CIFDocument* read_input(Input&& in) {
      CIFDocument* doc = new CIFDocument();
      parse_input(*doc, in);
      return doc;
    }
    
    inline CIFDocument* read_file(const std::string& filename) {
      pegtl::file_input<> in(filename);
      return read_input(in);
    }
    
    inline CIFDocument* read_string(const std::string& data) {
      pegtl::memory_input<> in(data, "string");
      return read_input(in);
    }
    
    inline CIFDocument* read_memory(const char* data, size_t size, const char* name) {
      pegtl::memory_input<> in(data, size, name);
      return read_input(in);
    }
    
    inline CIFDocument* read_cstream(std::FILE *f, size_t maximum, const char* name) {
      pegtl::cstream_input<> in(f, maximum, name);
      return read_input(in);
    }
    
    inline CIFDocument* read_istream(std::istream &is,
                                    size_t maximum,
                                    const char* name) {
      pegtl::istream_input<> in(is, maximum, name);
      return read_input(in);
    }
    */
    
    template<typename Input>
    bool check_syntax(Input&& in, std::string* msg) {
      try {
        return pegtl::parse<rules::file, pegtl::nothing, Errors>(in);
      } catch (pegtl::parse_error& e) {
        if (msg)
          *msg = e.what();
        return false;
      }
    }
    
    // A function for transparent reading of stdin and/or gzipped files
    // in addition to normal CIF files. For example, if used as:
    //   Document doc = read(MaybeStdin(argv[1]));
    // it reads from stdin if argv[1] is "-", or from a file otherwise.
    template<typename T>
    void read(T&& input, CIFDocument& d) {
      if (input.is_stdin())
        //return read_cstream(stdin, 16*1024, "stdin", d);
        read_cstream(stdin, 16*1024, "stdin", d);
      if (std::unique_ptr<char[]> mem = input.memory())
        //return read_memory(mem.get(), input.mem_size(), input.path().c_str());
        read_memory(mem.get(), input.mem_size(), input.path().c_str(), d);
      //return read_file(input.path());
      read_file(input.path(), d);
    }
    
    template<typename T>
    bool check_syntax_any(T&& input, std::string* msg) {
      if (std::unique_ptr<char[]> mem = input.memory()) {
        pegtl::memory_input<> in(mem.get(), input.mem_size(), input.path());
        return check_syntax(in, msg);
      }
      pegtl::file_input<> in(input.path());
      return check_syntax(in, msg);
    }
  } // namespace mmcif
} // namespace bilab
#endif /* mmcif_h */
