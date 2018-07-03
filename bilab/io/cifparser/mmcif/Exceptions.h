//
//  Exceptions.h
//  cifparse
//
//  Created by 曹巍 on 2018/06/21.
//  Copyright © 2018年 巍 曹. All rights reserved.
//

#ifndef Exceptions_h
#define Exceptions_h

#include <stdexcept>
#include <string>

namespace bilab {
  namespace mmcif {
    /// Base class for all RCSB exceptions
    class RcsbException : public std::exception
    {
    protected:
      std::string _message;
      
    public:
      RcsbException(const std::string& message = std::string(),
                    const std::string& location = std::string()){
        AppendMessage(message, location);
      }
      ~RcsbException() throw() {
        _message.clear();
      }
      
      void AppendMessage(const std::string& message = std::string(),
                         const std::string& location = std::string()){
        _message += "Message: \"" + message + "\" generated at: " +
          location + '\n';
      }
      
      const char* what() const throw(){
        return _message.c_str();
      }
    };
    
    
    /// Empty value exception (e.g. NULL pointer, empty string)
    class EmptyValueException : public RcsbException
    {
      
    public:
      EmptyValueException(const std::string& message = std::string(),
                          const std::string& location = std::string())
      {}
      
      ~EmptyValueException() throw() {};
      
    };
    
    
    /// Object not found (thrown everywhere except from .find() methods)
    class NotFoundException : public RcsbException
    {
      
    public:
      NotFoundException(const std::string& message = std::string(),
                        const std::string& location = std::string())
      {}
      
      ~NotFoundException() throw()
      {}

    };
    
    
    /// Object already exists
    class AlreadyExistsException : public RcsbException
    {
      
    public:
      AlreadyExistsException(const std::string& message = std::string(),
                             const std::string& location = std::string())
      {}
      ~AlreadyExistsException() throw()
      {}
      
    };
    
    /// Empty container
    class EmptyContainerException : public RcsbException
    {
      
    public:
      EmptyContainerException(const std::string& message = std::string(),
                              const std::string& location = std::string());
      ~EmptyContainerException() throw();
      
    };
    
    
    /// File mode exception (e.g. attempt to write to read-only file, invalid mode.)
    class FileModeException : public RcsbException
    {
      
    public:
      FileModeException(const std::string& message = std::string(),
                        const std::string& location = std::string())
      {}
      
      ~FileModeException() throw()
      {}
    };
    
    
    /// Invalid state exception (e.g. getting a row reference in a column-wise table/// )
    class InvalidStateException : public RcsbException
    {
      
    public:
      InvalidStateException(const std::string& message = std::string(),
                            const std::string& location = std::string())
      {}
      ~InvalidStateException() throw()
      {}
      
    };
    
    
    /// Generic files related exception (e.g. read error, write errror, etc.)
    class FileException : public RcsbException
    {
      
    public:
      FileException(const std::string& message = std::string(),
                    const std::string& location = std::string())
      {}
      ~FileException() throw()
      {}
      
    };
    
    /// Invalid command line options
    class InvalidOptionsException : public RcsbException
    {
    public:
      InvalidOptionsException(const std::string& message = std::string(),
                              const std::string& location = std::string())
      {
        _message = message;
      }
      ~InvalidOptionsException() throw()
      {}
    };
    
    /// Versions do not match
    class VersionMismatchException : public RcsbException
    {
      
    public:
      VersionMismatchException(const std::string& message = std::string(),
                               const std::string& location = std::string())
      {}
      ~VersionMismatchException() throw()
      {}
    };

  }
}
#endif /* Exceptions_h */
