//
//  utils.h
//  cifparse
//
//  Created by 曹巍 on 2018/06/21.
//  Copyright © 2018年 巍 曹. All rights reserved.
//

#ifndef utils_h
#define utils_h

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace bilab {
  namespace mmcif {
    // Useful functions
    inline void assert_tag(const std::string& tag);
    
    inline void ensure_mmcif_category(std::string& cat);
    
    inline bool is_null(const std::string& value);
    
    inline std::string as_string(const std::string& value);
    
    inline std::string as_string(const std::string* value);
    
    inline char as_char(const std::string& value, char null);
    
    inline bool is_text_field(const std::string& val);
    
    inline std::string quote(std::string v);
    
    // Case-insensitive version. Assumes the prefix/suffix is lowercase and ascii.
    inline bool istarts_with(const std::string& str, const std::string& prefix);
    
    inline bool iends_with(const std::string& str, const std::string& suffix);
    
    inline bool giends_with(const std::string& str, const std::string& suffix);
    
    inline std::string to_lower(std::string str);
    
    inline std::string trim_str(const std::string& str);
    
    [[noreturn]]
    inline void fail(const std::string& msg);
    
    //[[noreturn]]
    inline void message(const std::string& msg);
    /**
     * @brief starts_with
     * @function starts_with mmcif.h
     *
     * @author Wei Cao caotiger@gmail.com
     * @date   2018-06-19
     */
    inline bool starts_with(const std::string& s, const std::string& prefix);
    /**
     * @brief ends_with
     * @function ends_with mmcif.h
     *
     * @author Wei Cao caotiger@gmail.com
     * @date   2018-06-19
     */
    inline bool ends_with(const std::string& s, const std::string& suffix);

    std::vector<std::string>
    split(const std::string& s, const std::string& delim, const bool keep_empty = true);
    
    std::vector<std::vector<std::string>> groupby(std::vector<std::string>& v,
                                                  int num_groups);
  } // namespace mmcif
} // namespace bilab
#endif /* utils_h */
