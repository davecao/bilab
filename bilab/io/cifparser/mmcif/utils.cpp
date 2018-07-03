//
//  utils.h
//  cifparse
//
//  Created by 曹巍 on 2018/06/19.
//  Copyright © 2018年 巍 曹. All rights reserved.
//

#include "utils.h"

namespace bilab {
  namespace mmcif {
    // Useful functions
    inline void assert_tag(const std::string& tag) {
      if (tag[0] != '_')
        std::cerr << "Tag should start with '_', got: "
        << tag
        << std::endl;
    }
    
    inline void ensure_mmcif_category(std::string& cat) {
      if (cat[0] != '_')
        std::cerr << "Category should start with '_', got: " << cat << '\n';
      if (*(cat.end() - 1) != '.')
        cat += '.';
    }
    
    inline bool is_null(const std::string& value) {
      return value.size() == 1 && (value[0] == '?' || value[0] == '.');
    }
    
    inline std::string as_string(const std::string& value) {
      if (value.empty() || is_null(value))
        return "";
      if (value[0] == '"' || value[0] == '\'')
        return std::string(value.begin() + 1, value.end() - 1);
      if (value[0] == ';' && value.size() > 2 && *(value.end() - 2) == '\n') {
        bool crlf = *(value.end() - 3) == '\r';
        return std::string(value.begin() + 1, value.end() - (crlf ? 3 : 2));
      }
      return value;
    }
    
    inline std::string as_string(const std::string* value) {
      return value ? as_string(*value) : std::string();
    }
    
    inline char as_char(const std::string& value, char null) {
      if (is_null(value))
        return null;
      if (value.size() < 2)
        return value[0];
      const std::string s = as_string(value);
      if (s.size() < 2)
        return s[0];
      std::cerr << "Not a single character: " << value <<"\n";
      return null; // Original code maybe cause a little problem.
    }
    
    inline bool is_text_field(const std::string& val) {
      size_t len = val.size();
      return len > 3 && val[0] == ';' && (val[len-2] == '\n' || val[len-2] == '\r');
    }
    
    inline std::string quote(std::string v) {
      // strings with '(' happen to be quoted in mmCIF files, so we do the same
      if (v.find_first_of(" \t\r\n'\"()[]{}#") == std::string::npos &&
          v[0] != '_' && v[0] != '$' && v[0] != '\0' &&
          (v.size() > 1 || (v[0] != '.' && v[0] != '?')))
        return v;
      if (std::memchr(v.c_str(), '\n', v.size()))
        return ";" + v + "\n;";
      if (std::memchr(v.c_str(), '\'', v.size()) == nullptr)
        return "'" + v + "'";
      if (std::memchr(v.c_str(), '"', v.size()) == nullptr)
        return '"' + v + '"';
      return ";" + v + "\n;";
    }
    // Case-insensitive version. Assumes the prefix/suffix is lowercase and ascii.
    inline bool istarts_with(const std::string& str, const std::string& prefix) {
      return str.length() >= prefix.length() &&
      std::equal(std::begin(prefix), std::end(prefix), str.begin(),
                 [](char c1, char c2) { return c1 == std::tolower(c2); });
    }
    inline bool iends_with(const std::string& str, const std::string& suffix) {
      size_t sl = suffix.length();
      return str.length() >= sl &&
      std::equal(std::begin(suffix), std::end(suffix), str.end() - sl,
                 [](char c1, char c2) { return c1 == std::tolower(c2); });
    }
    
    inline bool giends_with(const std::string& str, const std::string& suffix) {
      return iends_with(str, suffix) || iends_with(str, suffix + ".gz");
    }
    
    inline std::string to_lower(std::string str) {
      for (char& c : str)
        if (c >= 'A' && c <= 'Z')
          c |= 0x20;
      return str;
    }
    
    inline std::string trim_str(const std::string& str)
    {
      std::string ws = " \r\n\t";
      std::string::size_type first = str.find_first_not_of(ws);
      if (first == std::string::npos)
        return std::string{};
      std::string::size_type last = str.find_last_not_of(ws);
      return str.substr(first, last - first + 1);
    }
    
    [[noreturn]]
    inline void fail(const std::string& msg) { throw std::runtime_error(msg); }
    
    //[[noreturn]]
    inline void message(const std::string& msg) {
      std::cout<<msg<<std::endl; }

    /**
     * @brief starts_with
     * @function starts_with mmcif.h
     *
     * @author Wei Cao caotiger@gmail.com
     * @date   2018-06-19
     */
    inline bool starts_with(const std::string& s, const std::string& prefix) {
      auto size = prefix.size();
      if (s.size() < size) return false;
      return std::equal(std::begin(prefix), std::end(prefix), std::begin(s));
    }
    
    /**
     * @brief ends_with
     * @function ends_with mmcif.h
     *
     * @author Wei Cao caotiger@gmail.com
     * @date   2018-06-19
     */
    inline bool ends_with(const std::string& s, const std::string& suffix)  {
      auto size = suffix.size();
      if (s.size() > size) return false;
      return std::equal(std::begin(suffix), std::end(suffix), std::begin(s));
    }

    std::vector<std::string>
    split(const std::string& s, const std::string& delim,
          const bool keep_empty) {
      std::vector<std::string> result;
      if (delim.empty()) {
        result.push_back(s);
        return result;
      }
      std::string::const_iterator substart = s.begin(), subend;
      while (true) {
        subend = search(substart, s.end(), delim.begin(), delim.end());
        std::string temp(substart, subend);
        if (keep_empty || !temp.empty()) {
          result.push_back(temp);
        }
        if (subend == s.end()) {
          break;
        }
        substart = subend + delim.size();
      }
      return result;
    }
    
    std::vector<std::vector<std::string>> groupby(std::vector<std::string>& v,
                                                  int num_groups) {
      std::vector<
        std::vector<std::string>
      > groups(size_t(std::ceil(v.size() / double(num_groups))));
      
      //Partition
      size_t vb = 0;
      size_t inc = 0;
      for (auto i = v.cbegin(); i != v.cend(); i += inc, ++vb){
        int d2e = int(std::distance(i, v.cend()));
        inc = std::min(d2e, num_groups);
        groups[vb].assign(i, i + inc);
      }
      return groups;
    }
  } // namespace mmcif
} // namespace bilab