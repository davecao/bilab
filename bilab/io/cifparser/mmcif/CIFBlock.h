//
//  CIFBlock.h
//  cifparse
//
//  Created by 曹巍 on 2018/06/21.
//  Copyright © 2018年 巍 曹. All rights reserved.
//

#ifndef CIFBlock_h
#define CIFBlock_h

#include "CIFCommon.h"
#include "CIFTable.h"

namespace bilab {
  namespace mmcif {
    /**
     * @brief CIFBlock
     * @class  CIFBlock mmcif.h
     *     Corresponding to "data_"
     *
     *  @author Wei Cao caotiger@gmail.com
     *  @date   2018-06-19
     */
    class CIFBlock
    {
    public:
      /*! Constructor */
      CIFBlock(){
        table = new CIFTable();
      };

      CIFBlock(const std::string& name_)
      : name(name_)
      {
        table = new CIFTable();
      }
      
      ~CIFBlock(){
        if (this->table != nullptr) {
          delete table;
        }
      };
      
      void add_table(std::vector<std::string>& tags,
                     std::vector<std::string>& values) {
        // add column names
        std::vector<std::string> items;
        for (auto c: tags) {
          std::string item = c.substr(c.find(".")+1);
          items.push_back(item);
        }
        table->AddColumn(items);
        
        // add values
        int num_items = static_cast<int>(items.size());
        std::vector<std::string>::iterator itr;
        std::vector<std::vector<std::string>> tuples;
        tuples = groupby(values, num_items);
        for (size_t i = 0; i < tuples.size(); i++){
          table->AddColumn(tuples[i]);
        }
      }
      
      void setName(std::string& n){
        name = n;
      }
      
      void show_table() {
        table->show_table();
      }
      // Fields
      std::string name;
      CIFTable* table;
    };
    
  }  // namespace mmcif
} // namespace bilab

#endif /* CIFBlock_h */
