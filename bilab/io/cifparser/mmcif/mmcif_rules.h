//
//  cif_rules.h
//  biTKCXX
//
//  Created by 曹巍 on 2018/06/20.
//  Copyright © 2018年 巍 曹. All rights reserved.
//

#ifndef cif_rules_h
#define cif_rules_h

#include "tao/pegtl.hpp"

namespace bilab {
  namespace pegtl = tao::pegtl;
  namespace mmcif {
    namespace rules {
      inline uint8_t lookup_table(char c) {
        static const uint8_t table[256] = {
          // 0  1  2  3  4  5  6  7  8  9  A  B  C  D  E  F
          0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 2, 0, 0, // 0
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 1
          2, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, // 2
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, // 3
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, // 4
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, // 5
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, // 6
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, // 7
          // 128-255
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
        };
        return table[static_cast<unsigned char>(c)];
      }
    
      template<int TableVal> struct lookup_char {
        using analyze_t = pegtl::analysis::generic<pegtl::analysis::rule_type::ANY>;
        template<typename Input> static bool match(Input& in) {
          if (!in.empty() && lookup_table(in.peek_char()) == TableVal) {
            if (TableVal == 2)  // this set includes new-line
              in.bump(1);
            else
              in.bump_in_this_line(1);
            return true;
          }
          return false;
        }
      };
    
      // (g) Character sets.
      // OrdinaryCharacter: ! % &  ()*+,-./0-9:  <=>?@A-Z[]  \ ^  `a-z{|}~
      using ordinary_char = lookup_char<1>;
    
      using ws_char = lookup_char<2>;
    
      // !"#$%&'()*+,-./0-9:;<=>?@A-Z[\]^_`a-z{|}~
      struct nonblank_ch : pegtl::range<'!', '~'> {};
    
      // ascii space is just before '!'
      struct anyprint_ch : pegtl::ranges<' ', '~', '\t'> {};
    
    
      // (f) White space and comments.
      struct comment : pegtl::if_must<pegtl::one<'#'>, pegtl::until<pegtl::eolf>> {};
      struct whitespace : pegtl::plus<pegtl::sor<ws_char, comment>> {};
      struct ws_or_eof : pegtl::sor<whitespace, pegtl::eof> {};

      // (b) Reserved words.
      struct str_data : TAOCPP_PEGTL_ISTRING("data_") {};
      struct str_loop : TAOCPP_PEGTL_ISTRING("loop_") {};
      struct str_global : TAOCPP_PEGTL_ISTRING("global_") {};
      struct str_save : TAOCPP_PEGTL_ISTRING("save_") {};
      struct str_stop : TAOCPP_PEGTL_ISTRING("stop_") {};
      struct keyword : pegtl::sor<str_data, str_loop, str_global, str_save, str_stop> {};
    
      // (e) Character strings and text fields.
      template<typename Q>
      struct endq: pegtl::seq<
            Q,
            pegtl::at<pegtl::sor<pegtl::one<' ','\n','\r','\t','#'>, pegtl::eof>>>
      {};
      
      // strict rule would be:
      // template <typename Q> struct quoted_tail : until<endq<Q>, anyprint_ch> {};
      // but it was relaxed after PDB accepted 5q1h with non-ascii character
      template<typename Q>
      struct quoted_tail : pegtl::until<endq<Q>, pegtl::not_one<'\n'>>
      {};
      
      template<typename Q>
      struct quoted : pegtl::if_must<Q, quoted_tail<Q>>
      {};
      
      struct singlequoted : quoted<pegtl::one<'\''>>
      {};
      
      struct doublequoted : quoted<pegtl::one<'"'>>
      {};
      
      struct field_sep : pegtl::seq<pegtl::bol, pegtl::one<';'>>
      {};
      
      // CIF 2.0 requires whitespace after text field, so it'd be:
      // until<endq<field_sep>> instead of until<field_sep>.
      struct textfield
      : pegtl::if_must<field_sep, pegtl::until<field_sep>>
      {};
      
      struct unquoted
      : pegtl::seq<
          pegtl::not_at<keyword>, pegtl::not_at<pegtl::one<'_','$','#'>>,
          pegtl::plus<nonblank_ch>>
      {};
    
      // (c) Tags and values.
    
      // (a) Basic structure of CIF.
      struct datablockname
        : pegtl::plus<nonblank_ch>
      {};
      
      struct datablockheading
        : pegtl::sor<
          pegtl::if_must<str_data, datablockname>,
          str_global>
      {};
      
      struct tag
      : pegtl::seq<pegtl::one<'_'>, pegtl::plus<nonblank_ch>>
      {};
      
      // unquoted value made of ordinary characters only - for a typical mmCIF file
      // it is faster to check it first even if we backtrack on some values_.
      struct simunq
        : pegtl::seq<pegtl::plus<ordinary_char>, pegtl::at<ws_char>>
      {};
      
      struct value
        : pegtl::sor<simunq, singlequoted, doublequoted, textfield, unquoted>
      {};
      
      struct loop_tag : tag {};
      
      struct loop_value : value {};
      
      struct loop_end
        : pegtl::opt<str_stop, ws_or_eof>
      {};
      
      struct loop
        : pegtl::if_must<
            str_loop, whitespace,
            pegtl::plus<pegtl::seq<loop_tag, whitespace, pegtl::discard>>,
            pegtl::sor<pegtl::plus<pegtl::seq<loop_value, ws_or_eof, pegtl::discard>>,
            // handle incorrect CIF with empty loop
            pegtl::at<pegtl::sor<str_loop, pegtl::eof>>>,
            loop_end>
      {};
      struct dataitem
        : pegtl::if_must<tag, whitespace, value, ws_or_eof, pegtl::discard>
      {};
      
      struct framename
        : pegtl::plus<nonblank_ch>
      {};
      
      struct endframe
        : str_save
      {};
      
      struct frame
        : pegtl::if_must<str_save, framename, whitespace,
              pegtl::star<pegtl::sor<dataitem, loop>>,
              endframe, ws_or_eof>
      {};
      
      struct datablock
        : pegtl::seq<datablockheading, ws_or_eof,
              pegtl::star<pegtl::sor<dataitem, loop, frame>>>
      {};
      
      struct file
        : pegtl::must<pegtl::opt<whitespace>, pegtl::star<datablock>, pegtl::eof>
      {};
    }  //  namespace rules
    
  }  //  namespace mmcif
  
  

} // namespace bilab

#endif /* cif_rules_h */
