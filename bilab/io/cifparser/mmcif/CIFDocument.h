//
//  CIFDocument.h
//  cifparse
//
//  Created by 曹巍 on 2018/06/21.
//  Copyright © 2018年 巍 曹. All rights reserved.
//

#ifndef CIFDocument_h
#define CIFDocument_h

#include "CIFCommon.h"
#include "CIFBlock.h"

namespace bilab {
  namespace mmcif {
    /**
     * @brief CIFDocument
     * @class  CIFDocument mmcif.h
     *
     *
     *  @author Wei Cao caotiger@gmail.com
     *  @date   2018-06-19
     */
    class CIFDocument
    {
    public:
      /*! Constructor */
      CIFDocument()
      : isloop(false), line_number(-1)
      {};

      ~CIFDocument(){
        if (!blocks.empty()) {
          for (unsigned int i = 0; i < blocks.size(); i++){
            delete blocks[i];
          }
        }
      };
      /*
      void clear() noexcept {
        source.clear();
        blocks.clear();
      }
       */
      CIFBlock* find_block(const std::string& name) {
        //for (CIFBlock* b : blocks)
        for(auto itr = blocks.begin(); itr != blocks.end(); ++itr ){
          if ((*itr)->name == name) return *itr;
        }
        return nullptr;
      }
      
      const CIFBlock* find_block(const std::string& name) const {
        return const_cast<CIFDocument*>(this)->find_block(name);
      }
      
      void add_block_loop(){
        if(!this->_tags.empty() && !this->_values.empty()){
          std::vector<std::string> results =
            bilab::mmcif::split(_tags[0].substr(1), ".");
          std::string block_name = results[0];
          CIFBlock* block = new CIFBlock();
          block->name = block_name;
          block->add_table(_tags, _values);
          blocks.emplace_back(block);
          
          // Clear _tags and _values
          this->_values.clear();
          this->_tags.clear();
        }
      }
      
      int getNumBlocks() {
        return blocks.size();
      }
      
      void show_block_names() {
        for(auto itr = blocks.begin(); itr != blocks.end(); ++itr ){
          std::cout << (*itr)->name <<std::endl;
        }
      }
      
      // public attributes
      std::vector<CIFBlock*> blocks;
      
      std::string source;
      std::string currTag;
      std::vector<std::string> _tags;
      std::vector<std::string> _values;
      bool isloop;
      int line_number;
    };
  }  // namespace mmcif
}  //namespace bilab
#endif /* CIFDocument_h */
