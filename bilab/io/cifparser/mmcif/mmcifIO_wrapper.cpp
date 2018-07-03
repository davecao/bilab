//
//  mmcifIO_wrapper.cpp
// 
//
//  Created by 曹巍 on 2018/06/19.
//  Copyright © 2018年 巍 曹. All rights reserved.
//
#include <stdio.h>
#include <iostream>

#include "mmcifIO_wrapper.h"

namespace mmcif = bilab::mmcif;

bool reader_mmcif(std::string& filename, mmcif::CIFDocument& doc) {
  return mmcif::read_file(filename, doc);
}