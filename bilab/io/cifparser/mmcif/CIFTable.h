//
//  CIFTable.h
//  cifparse
//
//  Created by 曹巍 on 2018/06/21.
//  Copyright © 2018年 巍 曹. All rights reserved.
//

#ifndef CIFTable_h
#define CIFTable_h

#include "CIFCommon.h"


namespace bilab {
  namespace mmcif {
    /**
     * @brief CIFTable
     * @class  CIFTable mmcif.h
     *
     *  A table of tuples
     *  Block: atom_site
     *  Columns:
     *     group_PDB,
     *     id
     *     type_symbol
     *     label_atom_id
     *     label_alt_id
     *     label_comp_id
     *     label_asym_id
     *     label_entity_id
     *     label_seq_id
     *     pdbx_PDB_ins_code
     *     Cartn_x
     *     Cartn_y
     *     Cartn_z
     *     occupancy
     *     B_iso_or_equiv
     *     pdbx_formal_charge
     *     auth_seq_id
     *     auth_comp_id
     *     auth_asym_id
     *     auth_atom_id
     *     pdbx_PDB_model_num
     *
     * ATOM   1    N  N   . SER A 1 1   ? -47.333 0.941   8.834   1.00 52.56  ? 71  SER P N   1
     *
     *
     *
     *  @author Wei Cao caotiger@gmail.com
     *  @date   2018-06-19
     */
    class CIFTable {
    public:
      /*! Constructor */
      CIFTable() : _numCols(0)
      {};
      
      CIFTable(const CIFTable& table){
        _numCols = table._numCols;
        
        // This statement only copies pointer addresses and not the tuple
        // vectors referenced by those pointers.
        _tuples = table._tuples;
        
        // Here, tuple vectors themselves are copied and previously set
        // pointers are properly initialized.
        
        for (unsigned int tupleI = 0; tupleI < _tuples.size(); ++tupleI) {
          _tuples[tupleI] = new std::vector<std::string>(*(table._tuples[tupleI]));
        }
      }
      /**
       *  Destructs a table.
       */
      ~CIFTable(){ Clear();};
      /**
       *  Deletes all the content from the table.
       *
       * Table has 0 columns and 0 tuples.
       */
      void Clear() {
        if (!_tuples.empty()) {
          for (unsigned int tupleI = 0; tupleI < _tuples.size(); tupleI++){
            _tuples[tupleI]->clear();
            delete _tuples[tupleI];
          }
          _tuples.clear();
        }
        _numCols = 0;
      }
      /**
       *  Copies a tuple table to another table (assignment operator).
       *
       * @parameter table  - reference to the source table
       *
       * @return Reference to the destination table
       */
      CIFTable& operator=(const CIFTable& table) {
        if (this != &table) {
          Clear();
          
          _numCols = table._numCols;
          
          // This statement only copies pointer addresses and not the tuple
          // vectors referenced by those pointers.
          _tuples = table._tuples;
          
          // Here, tuple vectors themselves are copied and previously set
          // pointers are properly initialized.
          
          for (unsigned int tupleI = 0; tupleI < _tuples.size();
               ++tupleI) {
            _tuples[tupleI] = new std::vector<std::string>(*(table._tuples[tupleI]));
          }
        }
        return(*this);
      }
      
      /**
       *  Retrieves the number of tuples in the table.
       *
       * @parameters: None
       *
       * @return The number of tuples in the table.
       *
       */
      inline unsigned int GetNumTuples() const {
        return static_cast<unsigned int>(_tuples.size());
      }
      
      inline unsigned int GetNumColumns() const{
        return _numCols;
      }

      unsigned int AddTuple(const std::vector<std::string>& tuple =
                            std::vector<std::string>()){
         InsertTuple(static_cast<unsigned int>(_tuples.size()), tuple);
         return static_cast<unsigned int>(_tuples.size());
       }
        
      void InsertTuple(const unsigned int tupleIndex,
                       const std::vector<std::string>& tuple = std::vector<std::string>()){
          InsertTuple(tupleIndex, tuple.begin(), tuple.end());
      }

      void InsertTuple(const unsigned int tupleIndex,
                       std::vector<std::string>::const_iterator tupleBeg,
                       std::vector<std::string>::const_iterator tupleEnd){
        if (tupleIndex > _tuples.size()){
          // Wrong tuple index.
          throw std::out_of_range("Invalid tuple index in CIFTable::InsertTuple");
        }
        
        if (tupleBeg > tupleEnd){
          // Wrong tuple range.
          throw std::out_of_range("Invalid tuple range in CIFTable::InsertTuple");
        }
        
        unsigned int tupleSize = static_cast<unsigned int>(tupleEnd - tupleBeg);
        
        if ((GetNumColumns() != 0) && (tupleSize > GetNumColumns())) {
          // Wrong tuple index.
          throw std::out_of_range("Invalid tuple size in CIFTable::InsertTuple");
        }
        
        if ((tupleSize != 0) && (tupleIndex != 0) &&
            (_tuples[tupleIndex - 1]->empty())) {
          // Previous tuple empty. Cannot insert, since rectangularity will be
          // violated.
          throw std::out_of_range("Previous tuple empty CIFTable::InsertTuple");
        }
        
        if (_numCols == 0) {
          // For empty table, new tuple defines the initial number of columns.
          _numCols = tupleSize;
        }
        
        _tuples.insert(_tuples.begin() + tupleIndex, new std::vector<std::string>);
        
        if (_tuples[tupleIndex]->empty()) {
          // If table tuple empty, insert empty cells.
          _tuples[tupleIndex]->insert(_tuples[tupleIndex]->begin(), _numCols,
                                      std::string());
        }
        
        std::copy(tupleBeg, tupleEnd, _tuples[tupleIndex]->begin());
      }
      
      void FillTuple(const unsigned int tupleIndex,
                     const std::vector<std::string>& tuple,
                     const unsigned int fromColIndex = 0) {
        if (tupleIndex >= _tuples.size()) {
          // Wrong tuple index.
          throw std::out_of_range("Invalid tuple index in CIFTable::FillTuple");
        }
        
        if (tuple.empty()){
          return;
        }
        
        if (fromColIndex > GetNumColumns()) {
          // Wrong from column index.
          throw std::out_of_range("Invalid from column index in CIFTable::FillTuple");
        }
        
        if ((GetNumColumns() != 0) && (tuple.size() > GetNumColumns() - fromColIndex)){
          // Wrong tuple index.
          throw std::out_of_range("Invalid tuple size in CIFTable::FillTuple");
        }
        
        std::copy(tuple.begin(), tuple.end(), _tuples[tupleIndex]->begin() +
             IntColIndex(fromColIndex));
      }
      
      void GetTuple(std::vector<std::string>& tuple,
                    const unsigned int tupleIndex,
                    const unsigned int fromColIndex, unsigned int toColIndex){
        tuple.clear();
        
        if (tupleIndex >= _tuples.size()){
          // Wrong tuple index.
          throw std::out_of_range("Invalid tuple index in CIFTable::GetTuple");
        }
        
        for (unsigned int colI = fromColIndex; colI < toColIndex; ++colI){
          tuple.push_back(operator()(tupleIndex, colI));
        }
      }
      
      const std::vector<std::string>& GetTuple(const unsigned int tupleIndex) {
        if (tupleIndex >= _tuples.size()){
          // Wrong tuple index.
          throw std::out_of_range("Invalid tuple index in CIFTable::GetTuple");
        }
        return(*(_tuples[tupleIndex]));
      }
      
      void ClearTuple(const unsigned int tupleIndex){
        if (tupleIndex >= _tuples.size()){
          // Wrong tuple index.
          throw std::out_of_range("Invalid tuple index in CIFTable::ClearTuple");
        }
        
        for (unsigned int colI = 0; colI < _tuples[tupleIndex]->size();
             ++colI) {
          (*_tuples[tupleIndex])[colI].clear();
        }
      }
      
      void DeleteTuple(const unsigned int tupleIndex){
        if (tupleIndex >= _tuples.size()){
          // Wrong tuple index.
          throw std::out_of_range("Invalid tuple index in CIFTable::DeleteTuple");
        }
        
        _tuples[tupleIndex]->clear();
        
        delete _tuples[tupleIndex];
        
        _tuples.erase(_tuples.begin() + tupleIndex);
        
        if (_tuples.empty()){
          _numCols = 0;
        }
      }
      
      unsigned int AddColumn(const std::vector<std::string>& col =
                             std::vector<std::string>()) {
        return(InsertColumn(GetNumColumns(), col));
      }
      
      unsigned int InsertColumn(const unsigned int atColIndex,
                                const std::vector<std::string>& col =
                                  std::vector<std::string>()) {
        InsertColumn(atColIndex, col.begin(), col.end());
        return(GetNumColumns());
      }
    
      void InsertColumn(const unsigned int atColIndex,
                      std::vector<std::string>::const_iterator colBeg,
                      std::vector<std::string>::const_iterator colEnd) {
        if (atColIndex > GetNumColumns()){
          // Wrong column index.
          throw std::out_of_range("Invalid column index in CIFTable::InsertColumn");
        }
      
        if (colBeg > colEnd){
          // Wrong tuple range.
          throw std::out_of_range("Invalid column range in CIFTable::InsertColumn");
        }
      
        unsigned int colSize = static_cast<unsigned int>(colEnd - colBeg);
        if (colSize == 0) {
          throw EmptyValueException("Empty column",
                                    "CIFTable::InsertColumn");
        }
      
        if ((!_tuples.empty()) && (colSize > _tuples.size())){
          throw std::out_of_range("Invalid column size in CIFTable::InsertColumn");
        }
      
        if ((GetNumColumns() == 0) && (colSize > _tuples.size())){
          unsigned int numTuples = static_cast<unsigned int>(colSize - _tuples.size());
          for (unsigned int tupleI = 0; tupleI < numTuples; ++tupleI){
            _tuples.push_back(new std::vector<std::string>);
          }
        }
        _numCols++;
        
        unsigned int intRowIndex = IntColIndex(atColIndex);
        
        for (unsigned int tupleI = 0; tupleI < _tuples.size(); ++tupleI) {
          _tuples[tupleI]->insert(_tuples[tupleI]->begin() + intRowIndex,
                                  1, std::string());
        }

        FillColumn(atColIndex, colBeg, colEnd);
      }
    
      void FillColumn(const unsigned int colIndex,
                      const std::vector<std::string>& col,
                      const unsigned int fromTupleIndex = 0){
          FillColumn(colIndex, col.begin(), col.end(), fromTupleIndex);
      }
      
      void FillColumn(const unsigned int colIndex,
                      std::vector<std::string>::const_iterator colBeg,
                      std::vector<std::string>::const_iterator colEnd,
                      const unsigned int fromTupleIndex = 0){
        if (colIndex >= GetNumColumns()){
          // Wrong column index.
          throw std::out_of_range("Invalid column index in CIFTable::FillColumn");
        }
        
        if (colBeg > colEnd){
          // Wrong tuple range.
          throw std::out_of_range("Invalid column range in CIFTable::FillColumn");
        }
        
        if (fromTupleIndex >= _tuples.size()){
          // Wrong tuple range.
          throw std::out_of_range("Invalid from tuple index in CIFTable::FillColumn");
        }
        
        unsigned int colSize = static_cast<unsigned int>(colEnd - colBeg);
        
        if (colSize == 0){
          return;
        }
        
        if (colSize > (_tuples.size() - fromTupleIndex)) {
          throw std::out_of_range("Invalid column size in CIFTable::FillColumn");
        }
        
        unsigned int intRowIndex = IntColIndex(colIndex);
        
        unsigned int tupleI = 0;
        for (std::vector<std::string>::const_iterator currIter = colBeg;
             currIter < colEnd;
             ++currIter, ++tupleI){
          (*_tuples[fromTupleIndex + tupleI])[intRowIndex] =
            *currIter;
        }
      }
      
      void GetColumn(std::vector<std::string>& col,
                     const unsigned int colIndex,
                     const unsigned int fromTupleIndex,
                     unsigned int toTupleIndex){
        col.clear();
        
        if (colIndex >= GetNumColumns()){
          // Wrong column index.
          throw std::out_of_range("Invalid column index in CIFTable::GetColumn");
        }
        
        if (fromTupleIndex >= toTupleIndex){
          // To index less than from index.
          throw std::out_of_range("Invalid tuple index range in CIFTable::GetColumn");
        }
        
        if (fromTupleIndex >= _tuples.size()){
          // Wrong from tuple index.
          throw std::out_of_range("Invalid from tuple index in CIFTable::GetColumn");
        }
        
        if (toTupleIndex > _tuples.size()){
          // Wrong to tuple index.
          throw std::out_of_range("Invalid to tuple index in CIFTable::GetColumn");
        }
        
        // toTupleIndex is inclusive, so one more spot for it needs to be reserved
        col.reserve(toTupleIndex - fromTupleIndex);
        
        unsigned int intRowIndex = IntColIndex(colIndex);
        
        for (unsigned int tupleI = fromTupleIndex; tupleI < toTupleIndex;
             tupleI++){
          col.push_back((*_tuples[tupleI])[intRowIndex]);
        }
      }
      
      void ClearColumn(const unsigned int colIndex) {
        if (colIndex >= GetNumColumns()){
          // Wrong column index.
          throw std::out_of_range("Invalid column index in CIFTable::ClearColumn");
        }
        
        unsigned int intRowIndex = IntColIndex(colIndex);
        
        for (unsigned int tupleI = 0; tupleI < _tuples.size(); ++tupleI)
        {
          (*_tuples[tupleI])[intRowIndex].clear();
        }

      }
    
      void DeleteColumn(const unsigned int colIndex){
        if (colIndex >= GetNumColumns()){
          // Wrong column index.
          throw std::out_of_range("Invalid column index in CIFTable::DeleteColumn");
        }
        _numCols--;
      }

      
      std::string& operator()(const unsigned int tupleIndex,
                              const unsigned int colIndex) {
        if (tupleIndex >= _tuples.size()){
          // Wrong tuple index.
          throw std::out_of_range("Invalid tuple index in CIFTable");
        }
        
        if (colIndex >= GetNumColumns()) {
          // Wrong column index.
          throw std::out_of_range("Invalid column index in CIFTable");
        }
        return((*_tuples[tupleIndex])[IntColIndex(colIndex)]);
      }
      
      const std::string& operator()(const unsigned int tupleIndex,
                                    const unsigned int colIndex) const {
        if (tupleIndex >= _tuples.size()){
          // Wrong tuple index.
          throw std::out_of_range("Invalid tuple index in CIFTable::GetCell");
        }
        
        if (colIndex >= GetNumColumns()){
          // Wrong column index.
          throw std::out_of_range("Invalid column index in CIFTable::GetCell");
        }
        return((*_tuples[tupleIndex])[IntColIndex(colIndex)]);
      }
      
      
      void show_table(){
        int rows = GetNumTuples();
        int cols = GetNumColumns();
        std::cout << "There are " << rows <<"x"<< cols <<" in the doc" << std::endl;
        for (unsigned int i = 0; i < rows; i++) {
          for (unsigned int j = 0; j < cols; j++) {
            std::cout<< (*_tuples[i])[IntColIndex(j)] <<",";
          }
          std::cout << std::endl;
        }
      }
    private:

      unsigned int _numCols;
      //std::vector<std::string> _column_names;
      std::vector<std::vector<std::string>*> _tuples;
      
      inline unsigned int IntColIndex(const unsigned int colIndex) const {
        return colIndex;
      }
    };
  }  // namespace mmcif
} // namespace bilab
#endif /* CIFTable_h */
