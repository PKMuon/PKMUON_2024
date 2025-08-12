// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME Object
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "/home/Aser_fxr_9999/PKMUON_2024/include/Object.hh"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_Track(void *p = nullptr);
   static void *newArray_Track(Long_t size, void *p);
   static void delete_Track(void *p);
   static void deleteArray_Track(void *p);
   static void destruct_Track(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Track*)
   {
      ::Track *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Track >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("Track", ::Track::Class_Version(), "", 43,
                  typeid(::Track), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Track::Dictionary, isa_proxy, 4,
                  sizeof(::Track) );
      instance.SetNew(&new_Track);
      instance.SetNewArray(&newArray_Track);
      instance.SetDelete(&delete_Track);
      instance.SetDeleteArray(&deleteArray_Track);
      instance.SetDestructor(&destruct_Track);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Track*)
   {
      return GenerateInitInstanceLocal(static_cast<::Track*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Track*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_Params(void *p = nullptr);
   static void *newArray_Params(Long_t size, void *p);
   static void delete_Params(void *p);
   static void deleteArray_Params(void *p);
   static void destruct_Params(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Params*)
   {
      ::Params *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Params >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("Params", ::Params::Class_Version(), "", 62,
                  typeid(::Params), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Params::Dictionary, isa_proxy, 4,
                  sizeof(::Params) );
      instance.SetNew(&new_Params);
      instance.SetNewArray(&newArray_Params);
      instance.SetDelete(&delete_Params);
      instance.SetDeleteArray(&deleteArray_Params);
      instance.SetDestructor(&destruct_Params);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Params*)
   {
      return GenerateInitInstanceLocal(static_cast<::Params*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Params*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_Edep(void *p = nullptr);
   static void *newArray_Edep(Long_t size, void *p);
   static void delete_Edep(void *p);
   static void deleteArray_Edep(void *p);
   static void destruct_Edep(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Edep*)
   {
      ::Edep *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Edep >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("Edep", ::Edep::Class_Version(), "", 81,
                  typeid(::Edep), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Edep::Dictionary, isa_proxy, 4,
                  sizeof(::Edep) );
      instance.SetNew(&new_Edep);
      instance.SetNewArray(&newArray_Edep);
      instance.SetDelete(&delete_Edep);
      instance.SetDeleteArray(&deleteArray_Edep);
      instance.SetDestructor(&destruct_Edep);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Edep*)
   {
      return GenerateInitInstanceLocal(static_cast<::Edep*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Edep*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_Process(void *p = nullptr);
   static void *newArray_Process(Long_t size, void *p);
   static void delete_Process(void *p);
   static void deleteArray_Process(void *p);
   static void destruct_Process(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Process*)
   {
      ::Process *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Process >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("Process", ::Process::Class_Version(), "", 101,
                  typeid(::Process), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Process::Dictionary, isa_proxy, 4,
                  sizeof(::Process) );
      instance.SetNew(&new_Process);
      instance.SetNewArray(&newArray_Process);
      instance.SetDelete(&delete_Process);
      instance.SetDeleteArray(&deleteArray_Process);
      instance.SetDestructor(&destruct_Process);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Process*)
   {
      return GenerateInitInstanceLocal(static_cast<::Process*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Process*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_Event(void *p = nullptr);
   static void *newArray_Event(Long_t size, void *p);
   static void delete_Event(void *p);
   static void deleteArray_Event(void *p);
   static void destruct_Event(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Event*)
   {
      ::Event *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Event >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("Event", ::Event::Class_Version(), "", 115,
                  typeid(::Event), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Event::Dictionary, isa_proxy, 4,
                  sizeof(::Event) );
      instance.SetNew(&new_Event);
      instance.SetNewArray(&newArray_Event);
      instance.SetDelete(&delete_Event);
      instance.SetDeleteArray(&deleteArray_Event);
      instance.SetDestructor(&destruct_Event);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Event*)
   {
      return GenerateInitInstanceLocal(static_cast<::Event*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Event*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr Track::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *Track::Class_Name()
{
   return "Track";
}

//______________________________________________________________________________
const char *Track::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Track*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int Track::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Track*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Track::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Track*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Track::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Track*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr Params::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *Params::Class_Name()
{
   return "Params";
}

//______________________________________________________________________________
const char *Params::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Params*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int Params::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Params*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Params::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Params*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Params::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Params*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr Edep::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *Edep::Class_Name()
{
   return "Edep";
}

//______________________________________________________________________________
const char *Edep::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Edep*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int Edep::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Edep*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Edep::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Edep*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Edep::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Edep*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr Process::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *Process::Class_Name()
{
   return "Process";
}

//______________________________________________________________________________
const char *Process::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Process*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int Process::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Process*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Process::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Process*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Process::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Process*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr Event::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *Event::Class_Name()
{
   return "Event";
}

//______________________________________________________________________________
const char *Event::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Event*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int Event::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Event*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Event::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Event*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Event::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Event*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void Track::Streamer(TBuffer &R__b)
{
   // Stream an object of class Track.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Track::Class(),this);
   } else {
      R__b.WriteClassBuffer(Track::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Track(void *p) {
      return  p ? new(p) ::Track : new ::Track;
   }
   static void *newArray_Track(Long_t nElements, void *p) {
      return p ? new(p) ::Track[nElements] : new ::Track[nElements];
   }
   // Wrapper around operator delete
   static void delete_Track(void *p) {
      delete (static_cast<::Track*>(p));
   }
   static void deleteArray_Track(void *p) {
      delete [] (static_cast<::Track*>(p));
   }
   static void destruct_Track(void *p) {
      typedef ::Track current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Track

//______________________________________________________________________________
void Params::Streamer(TBuffer &R__b)
{
   // Stream an object of class Params.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Params::Class(),this);
   } else {
      R__b.WriteClassBuffer(Params::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Params(void *p) {
      return  p ? new(p) ::Params : new ::Params;
   }
   static void *newArray_Params(Long_t nElements, void *p) {
      return p ? new(p) ::Params[nElements] : new ::Params[nElements];
   }
   // Wrapper around operator delete
   static void delete_Params(void *p) {
      delete (static_cast<::Params*>(p));
   }
   static void deleteArray_Params(void *p) {
      delete [] (static_cast<::Params*>(p));
   }
   static void destruct_Params(void *p) {
      typedef ::Params current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Params

//______________________________________________________________________________
void Edep::Streamer(TBuffer &R__b)
{
   // Stream an object of class Edep.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Edep::Class(),this);
   } else {
      R__b.WriteClassBuffer(Edep::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Edep(void *p) {
      return  p ? new(p) ::Edep : new ::Edep;
   }
   static void *newArray_Edep(Long_t nElements, void *p) {
      return p ? new(p) ::Edep[nElements] : new ::Edep[nElements];
   }
   // Wrapper around operator delete
   static void delete_Edep(void *p) {
      delete (static_cast<::Edep*>(p));
   }
   static void deleteArray_Edep(void *p) {
      delete [] (static_cast<::Edep*>(p));
   }
   static void destruct_Edep(void *p) {
      typedef ::Edep current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Edep

//______________________________________________________________________________
void Process::Streamer(TBuffer &R__b)
{
   // Stream an object of class Process.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Process::Class(),this);
   } else {
      R__b.WriteClassBuffer(Process::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Process(void *p) {
      return  p ? new(p) ::Process : new ::Process;
   }
   static void *newArray_Process(Long_t nElements, void *p) {
      return p ? new(p) ::Process[nElements] : new ::Process[nElements];
   }
   // Wrapper around operator delete
   static void delete_Process(void *p) {
      delete (static_cast<::Process*>(p));
   }
   static void deleteArray_Process(void *p) {
      delete [] (static_cast<::Process*>(p));
   }
   static void destruct_Process(void *p) {
      typedef ::Process current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Process

//______________________________________________________________________________
void Event::Streamer(TBuffer &R__b)
{
   // Stream an object of class Event.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Event::Class(),this);
   } else {
      R__b.WriteClassBuffer(Event::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Event(void *p) {
      return  p ? new(p) ::Event : new ::Event;
   }
   static void *newArray_Event(Long_t nElements, void *p) {
      return p ? new(p) ::Event[nElements] : new ::Event[nElements];
   }
   // Wrapper around operator delete
   static void delete_Event(void *p) {
      delete (static_cast<::Event*>(p));
   }
   static void deleteArray_Event(void *p) {
      delete [] (static_cast<::Event*>(p));
   }
   static void destruct_Event(void *p) {
      typedef ::Event current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Event

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = nullptr);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 389,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<double>","std::vector<double, std::allocator<double> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr))->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete (static_cast<vector<double>*>(p));
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] (static_cast<vector<double>*>(p));
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   // Registration Schema evolution read functions
   int RecordReadRules_Object() {
      return 0;
   }
   static int _R__UNIQUE_DICT_(ReadRules_Object) = RecordReadRules_Object();R__UseDummy(_R__UNIQUE_DICT_(ReadRules_Object));
} // namespace ROOT
namespace {
  void TriggerDictionaryInitialization_Object_Impl() {
    static const char* headers[] = {
"/home/Aser_fxr_9999/PKMUON_2024/include/Object.hh",
nullptr
    };
    static const char* includePaths[] = {
"/home/Aser_fxr_9999/root/include/",
"/home/Aser_fxr_9999/PKMUON_2024/build1/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "Object dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$/home/Aser_fxr_9999/PKMUON_2024/include/Object.hh")))  Track;
class __attribute__((annotate("$clingAutoload$/home/Aser_fxr_9999/PKMUON_2024/include/Object.hh")))  Params;
class __attribute__((annotate("$clingAutoload$/home/Aser_fxr_9999/PKMUON_2024/include/Object.hh")))  Edep;
class __attribute__((annotate("$clingAutoload$/home/Aser_fxr_9999/PKMUON_2024/include/Object.hh")))  Process;
class __attribute__((annotate("$clingAutoload$/home/Aser_fxr_9999/PKMUON_2024/include/Object.hh")))  Event;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "Object dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "/home/Aser_fxr_9999/PKMUON_2024/include/Object.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"Edep", payloadCode, "@",
"Event", payloadCode, "@",
"Params", payloadCode, "@",
"Process", payloadCode, "@",
"Track", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("Object",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_Object_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_Object_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_Object() {
  TriggerDictionaryInitialization_Object_Impl();
}
