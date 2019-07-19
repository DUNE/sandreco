// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME struct

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
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

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "struct.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *cell_Dictionary();
   static void cell_TClassManip(TClass*);
   static void *new_cell(void *p = 0);
   static void *newArray_cell(Long_t size, void *p);
   static void delete_cell(void *p);
   static void deleteArray_cell(void *p);
   static void destruct_cell(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::cell*)
   {
      ::cell *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::cell));
      static ::ROOT::TGenericClassInfo 
         instance("cell", "struct.h", 16,
                  typeid(::cell), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &cell_Dictionary, isa_proxy, 4,
                  sizeof(::cell) );
      instance.SetNew(&new_cell);
      instance.SetNewArray(&newArray_cell);
      instance.SetDelete(&delete_cell);
      instance.SetDeleteArray(&deleteArray_cell);
      instance.SetDestructor(&destruct_cell);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::cell*)
   {
      return GenerateInitInstanceLocal((::cell*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::cell*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *cell_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::cell*)0x0)->GetClass();
      cell_TClassManip(theClass);
   return theClass;
   }

   static void cell_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *cluster_Dictionary();
   static void cluster_TClassManip(TClass*);
   static void *new_cluster(void *p = 0);
   static void *newArray_cluster(Long_t size, void *p);
   static void delete_cluster(void *p);
   static void deleteArray_cluster(void *p);
   static void destruct_cluster(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::cluster*)
   {
      ::cluster *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::cluster));
      static ::ROOT::TGenericClassInfo 
         instance("cluster", "struct.h", 35,
                  typeid(::cluster), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &cluster_Dictionary, isa_proxy, 4,
                  sizeof(::cluster) );
      instance.SetNew(&new_cluster);
      instance.SetNewArray(&newArray_cluster);
      instance.SetDelete(&delete_cluster);
      instance.SetDeleteArray(&deleteArray_cluster);
      instance.SetDestructor(&destruct_cluster);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::cluster*)
   {
      return GenerateInitInstanceLocal((::cluster*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::cluster*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *cluster_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::cluster*)0x0)->GetClass();
      cluster_TClassManip(theClass);
   return theClass;
   }

   static void cluster_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *digit_Dictionary();
   static void digit_TClassManip(TClass*);
   static void *new_digit(void *p = 0);
   static void *newArray_digit(Long_t size, void *p);
   static void delete_digit(void *p);
   static void deleteArray_digit(void *p);
   static void destruct_digit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::digit*)
   {
      ::digit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::digit));
      static ::ROOT::TGenericClassInfo 
         instance("digit", "struct.h", 66,
                  typeid(::digit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &digit_Dictionary, isa_proxy, 4,
                  sizeof(::digit) );
      instance.SetNew(&new_digit);
      instance.SetNewArray(&newArray_digit);
      instance.SetDelete(&delete_digit);
      instance.SetDeleteArray(&deleteArray_digit);
      instance.SetDestructor(&destruct_digit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::digit*)
   {
      return GenerateInitInstanceLocal((::digit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::digit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *digit_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::digit*)0x0)->GetClass();
      digit_TClassManip(theClass);
   return theClass;
   }

   static void digit_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *track_Dictionary();
   static void track_TClassManip(TClass*);
   static void *new_track(void *p = 0);
   static void *newArray_track(Long_t size, void *p);
   static void delete_track(void *p);
   static void deleteArray_track(void *p);
   static void destruct_track(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::track*)
   {
      ::track *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::track));
      static ::ROOT::TGenericClassInfo 
         instance("track", "struct.h", 77,
                  typeid(::track), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &track_Dictionary, isa_proxy, 4,
                  sizeof(::track) );
      instance.SetNew(&new_track);
      instance.SetNewArray(&newArray_track);
      instance.SetDelete(&delete_track);
      instance.SetDeleteArray(&deleteArray_track);
      instance.SetDestructor(&destruct_track);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::track*)
   {
      return GenerateInitInstanceLocal((::track*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::track*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *track_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::track*)0x0)->GetClass();
      track_TClassManip(theClass);
   return theClass;
   }

   static void track_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *particle_Dictionary();
   static void particle_TClassManip(TClass*);
   static void *new_particle(void *p = 0);
   static void *newArray_particle(Long_t size, void *p);
   static void delete_particle(void *p);
   static void deleteArray_particle(void *p);
   static void destruct_particle(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::particle*)
   {
      ::particle *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::particle));
      static ::ROOT::TGenericClassInfo 
         instance("particle", "struct.h", 96,
                  typeid(::particle), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &particle_Dictionary, isa_proxy, 4,
                  sizeof(::particle) );
      instance.SetNew(&new_particle);
      instance.SetNewArray(&newArray_particle);
      instance.SetDelete(&delete_particle);
      instance.SetDeleteArray(&deleteArray_particle);
      instance.SetDestructor(&destruct_particle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::particle*)
   {
      return GenerateInitInstanceLocal((::particle*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::particle*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *particle_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::particle*)0x0)->GetClass();
      particle_TClassManip(theClass);
   return theClass;
   }

   static void particle_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *event_Dictionary();
   static void event_TClassManip(TClass*);
   static void *new_event(void *p = 0);
   static void *newArray_event(Long_t size, void *p);
   static void delete_event(void *p);
   static void deleteArray_event(void *p);
   static void destruct_event(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::event*)
   {
      ::event *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::event));
      static ::ROOT::TGenericClassInfo 
         instance("event", "struct.h", 132,
                  typeid(::event), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &event_Dictionary, isa_proxy, 4,
                  sizeof(::event) );
      instance.SetNew(&new_event);
      instance.SetNewArray(&newArray_event);
      instance.SetDelete(&delete_event);
      instance.SetDeleteArray(&deleteArray_event);
      instance.SetDestructor(&destruct_event);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::event*)
   {
      return GenerateInitInstanceLocal((::event*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::event*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *event_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::event*)0x0)->GetClass();
      event_TClassManip(theClass);
   return theClass;
   }

   static void event_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_cell(void *p) {
      return  p ? new(p) ::cell : new ::cell;
   }
   static void *newArray_cell(Long_t nElements, void *p) {
      return p ? new(p) ::cell[nElements] : new ::cell[nElements];
   }
   // Wrapper around operator delete
   static void delete_cell(void *p) {
      delete ((::cell*)p);
   }
   static void deleteArray_cell(void *p) {
      delete [] ((::cell*)p);
   }
   static void destruct_cell(void *p) {
      typedef ::cell current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::cell

namespace ROOT {
   // Wrappers around operator new
   static void *new_cluster(void *p) {
      return  p ? new(p) ::cluster : new ::cluster;
   }
   static void *newArray_cluster(Long_t nElements, void *p) {
      return p ? new(p) ::cluster[nElements] : new ::cluster[nElements];
   }
   // Wrapper around operator delete
   static void delete_cluster(void *p) {
      delete ((::cluster*)p);
   }
   static void deleteArray_cluster(void *p) {
      delete [] ((::cluster*)p);
   }
   static void destruct_cluster(void *p) {
      typedef ::cluster current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::cluster

namespace ROOT {
   // Wrappers around operator new
   static void *new_digit(void *p) {
      return  p ? new(p) ::digit : new ::digit;
   }
   static void *newArray_digit(Long_t nElements, void *p) {
      return p ? new(p) ::digit[nElements] : new ::digit[nElements];
   }
   // Wrapper around operator delete
   static void delete_digit(void *p) {
      delete ((::digit*)p);
   }
   static void deleteArray_digit(void *p) {
      delete [] ((::digit*)p);
   }
   static void destruct_digit(void *p) {
      typedef ::digit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::digit

namespace ROOT {
   // Wrappers around operator new
   static void *new_track(void *p) {
      return  p ? new(p) ::track : new ::track;
   }
   static void *newArray_track(Long_t nElements, void *p) {
      return p ? new(p) ::track[nElements] : new ::track[nElements];
   }
   // Wrapper around operator delete
   static void delete_track(void *p) {
      delete ((::track*)p);
   }
   static void deleteArray_track(void *p) {
      delete [] ((::track*)p);
   }
   static void destruct_track(void *p) {
      typedef ::track current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::track

namespace ROOT {
   // Wrappers around operator new
   static void *new_particle(void *p) {
      return  p ? new(p) ::particle : new ::particle;
   }
   static void *newArray_particle(Long_t nElements, void *p) {
      return p ? new(p) ::particle[nElements] : new ::particle[nElements];
   }
   // Wrapper around operator delete
   static void delete_particle(void *p) {
      delete ((::particle*)p);
   }
   static void deleteArray_particle(void *p) {
      delete [] ((::particle*)p);
   }
   static void destruct_particle(void *p) {
      typedef ::particle current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::particle

namespace ROOT {
   // Wrappers around operator new
   static void *new_event(void *p) {
      return  p ? new(p) ::event : new ::event;
   }
   static void *newArray_event(Long_t nElements, void *p) {
      return p ? new(p) ::event[nElements] : new ::event[nElements];
   }
   // Wrapper around operator delete
   static void delete_event(void *p) {
      delete ((::event*)p);
   }
   static void deleteArray_event(void *p) {
      delete [] ((::event*)p);
   }
   static void destruct_event(void *p) {
      typedef ::event current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::event

namespace ROOT {
   static TClass *vectorlEtrackgR_Dictionary();
   static void vectorlEtrackgR_TClassManip(TClass*);
   static void *new_vectorlEtrackgR(void *p = 0);
   static void *newArray_vectorlEtrackgR(Long_t size, void *p);
   static void delete_vectorlEtrackgR(void *p);
   static void deleteArray_vectorlEtrackgR(void *p);
   static void destruct_vectorlEtrackgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<track>*)
   {
      vector<track> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<track>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<track>", -2, "vector", 210,
                  typeid(vector<track>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEtrackgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<track>) );
      instance.SetNew(&new_vectorlEtrackgR);
      instance.SetNewArray(&newArray_vectorlEtrackgR);
      instance.SetDelete(&delete_vectorlEtrackgR);
      instance.SetDeleteArray(&deleteArray_vectorlEtrackgR);
      instance.SetDestructor(&destruct_vectorlEtrackgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<track> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<track>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEtrackgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<track>*)0x0)->GetClass();
      vectorlEtrackgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEtrackgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEtrackgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<track> : new vector<track>;
   }
   static void *newArray_vectorlEtrackgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<track>[nElements] : new vector<track>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEtrackgR(void *p) {
      delete ((vector<track>*)p);
   }
   static void deleteArray_vectorlEtrackgR(void *p) {
      delete [] ((vector<track>*)p);
   }
   static void destruct_vectorlEtrackgR(void *p) {
      typedef vector<track> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<track>

namespace ROOT {
   static TClass *vectorlEparticlegR_Dictionary();
   static void vectorlEparticlegR_TClassManip(TClass*);
   static void *new_vectorlEparticlegR(void *p = 0);
   static void *newArray_vectorlEparticlegR(Long_t size, void *p);
   static void delete_vectorlEparticlegR(void *p);
   static void deleteArray_vectorlEparticlegR(void *p);
   static void destruct_vectorlEparticlegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<particle>*)
   {
      vector<particle> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<particle>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<particle>", -2, "vector", 210,
                  typeid(vector<particle>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEparticlegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<particle>) );
      instance.SetNew(&new_vectorlEparticlegR);
      instance.SetNewArray(&newArray_vectorlEparticlegR);
      instance.SetDelete(&delete_vectorlEparticlegR);
      instance.SetDeleteArray(&deleteArray_vectorlEparticlegR);
      instance.SetDestructor(&destruct_vectorlEparticlegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<particle> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<particle>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEparticlegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<particle>*)0x0)->GetClass();
      vectorlEparticlegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEparticlegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEparticlegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<particle> : new vector<particle>;
   }
   static void *newArray_vectorlEparticlegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<particle>[nElements] : new vector<particle>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEparticlegR(void *p) {
      delete ((vector<particle>*)p);
   }
   static void deleteArray_vectorlEparticlegR(void *p) {
      delete [] ((vector<particle>*)p);
   }
   static void destruct_vectorlEparticlegR(void *p) {
      typedef vector<particle> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<particle>

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 210,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEhitgR_Dictionary();
   static void vectorlEhitgR_TClassManip(TClass*);
   static void *new_vectorlEhitgR(void *p = 0);
   static void *newArray_vectorlEhitgR(Long_t size, void *p);
   static void delete_vectorlEhitgR(void *p);
   static void deleteArray_vectorlEhitgR(void *p);
   static void destruct_vectorlEhitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<hit>*)
   {
      vector<hit> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<hit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<hit>", -2, "vector", 210,
                  typeid(vector<hit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEhitgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<hit>) );
      instance.SetNew(&new_vectorlEhitgR);
      instance.SetNewArray(&newArray_vectorlEhitgR);
      instance.SetDelete(&delete_vectorlEhitgR);
      instance.SetDeleteArray(&deleteArray_vectorlEhitgR);
      instance.SetDestructor(&destruct_vectorlEhitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<hit> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<hit>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEhitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<hit>*)0x0)->GetClass();
      vectorlEhitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEhitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEhitgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<hit> : new vector<hit>;
   }
   static void *newArray_vectorlEhitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<hit>[nElements] : new vector<hit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEhitgR(void *p) {
      delete ((vector<hit>*)p);
   }
   static void deleteArray_vectorlEhitgR(void *p) {
      delete [] ((vector<hit>*)p);
   }
   static void destruct_vectorlEhitgR(void *p) {
      typedef vector<hit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<hit>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 210,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlEdigitgR_Dictionary();
   static void vectorlEdigitgR_TClassManip(TClass*);
   static void *new_vectorlEdigitgR(void *p = 0);
   static void *newArray_vectorlEdigitgR(Long_t size, void *p);
   static void delete_vectorlEdigitgR(void *p);
   static void deleteArray_vectorlEdigitgR(void *p);
   static void destruct_vectorlEdigitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<digit>*)
   {
      vector<digit> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<digit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<digit>", -2, "vector", 210,
                  typeid(vector<digit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdigitgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<digit>) );
      instance.SetNew(&new_vectorlEdigitgR);
      instance.SetNewArray(&newArray_vectorlEdigitgR);
      instance.SetDelete(&delete_vectorlEdigitgR);
      instance.SetDeleteArray(&deleteArray_vectorlEdigitgR);
      instance.SetDestructor(&destruct_vectorlEdigitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<digit> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<digit>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdigitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<digit>*)0x0)->GetClass();
      vectorlEdigitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdigitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdigitgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<digit> : new vector<digit>;
   }
   static void *newArray_vectorlEdigitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<digit>[nElements] : new vector<digit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdigitgR(void *p) {
      delete ((vector<digit>*)p);
   }
   static void deleteArray_vectorlEdigitgR(void *p) {
      delete [] ((vector<digit>*)p);
   }
   static void destruct_vectorlEdigitgR(void *p) {
      typedef vector<digit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<digit>

namespace ROOT {
   static TClass *vectorlEclustergR_Dictionary();
   static void vectorlEclustergR_TClassManip(TClass*);
   static void *new_vectorlEclustergR(void *p = 0);
   static void *newArray_vectorlEclustergR(Long_t size, void *p);
   static void delete_vectorlEclustergR(void *p);
   static void deleteArray_vectorlEclustergR(void *p);
   static void destruct_vectorlEclustergR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<cluster>*)
   {
      vector<cluster> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<cluster>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<cluster>", -2, "vector", 210,
                  typeid(vector<cluster>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEclustergR_Dictionary, isa_proxy, 4,
                  sizeof(vector<cluster>) );
      instance.SetNew(&new_vectorlEclustergR);
      instance.SetNewArray(&newArray_vectorlEclustergR);
      instance.SetDelete(&delete_vectorlEclustergR);
      instance.SetDeleteArray(&deleteArray_vectorlEclustergR);
      instance.SetDestructor(&destruct_vectorlEclustergR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<cluster> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<cluster>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEclustergR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<cluster>*)0x0)->GetClass();
      vectorlEclustergR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEclustergR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEclustergR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<cluster> : new vector<cluster>;
   }
   static void *newArray_vectorlEclustergR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<cluster>[nElements] : new vector<cluster>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEclustergR(void *p) {
      delete ((vector<cluster>*)p);
   }
   static void deleteArray_vectorlEclustergR(void *p) {
      delete [] ((vector<cluster>*)p);
   }
   static void destruct_vectorlEclustergR(void *p) {
      typedef vector<cluster> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<cluster>

namespace ROOT {
   static TClass *vectorlEcellgR_Dictionary();
   static void vectorlEcellgR_TClassManip(TClass*);
   static void *new_vectorlEcellgR(void *p = 0);
   static void *newArray_vectorlEcellgR(Long_t size, void *p);
   static void delete_vectorlEcellgR(void *p);
   static void deleteArray_vectorlEcellgR(void *p);
   static void destruct_vectorlEcellgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<cell>*)
   {
      vector<cell> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<cell>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<cell>", -2, "vector", 210,
                  typeid(vector<cell>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEcellgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<cell>) );
      instance.SetNew(&new_vectorlEcellgR);
      instance.SetNewArray(&newArray_vectorlEcellgR);
      instance.SetDelete(&delete_vectorlEcellgR);
      instance.SetDeleteArray(&deleteArray_vectorlEcellgR);
      instance.SetDestructor(&destruct_vectorlEcellgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<cell> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<cell>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEcellgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<cell>*)0x0)->GetClass();
      vectorlEcellgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEcellgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEcellgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<cell> : new vector<cell>;
   }
   static void *newArray_vectorlEcellgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<cell>[nElements] : new vector<cell>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEcellgR(void *p) {
      delete ((vector<cell>*)p);
   }
   static void deleteArray_vectorlEcellgR(void *p) {
      delete [] ((vector<cell>*)p);
   }
   static void destruct_vectorlEcellgR(void *p) {
      typedef vector<cell> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<cell>

namespace ROOT {
   static TClass *maplEstringcOvectorlEhitgRsPgR_Dictionary();
   static void maplEstringcOvectorlEhitgRsPgR_TClassManip(TClass*);
   static void *new_maplEstringcOvectorlEhitgRsPgR(void *p = 0);
   static void *newArray_maplEstringcOvectorlEhitgRsPgR(Long_t size, void *p);
   static void delete_maplEstringcOvectorlEhitgRsPgR(void *p);
   static void deleteArray_maplEstringcOvectorlEhitgRsPgR(void *p);
   static void destruct_maplEstringcOvectorlEhitgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,vector<hit> >*)
   {
      map<string,vector<hit> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,vector<hit> >));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,vector<hit> >", -2, "map", 96,
                  typeid(map<string,vector<hit> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOvectorlEhitgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(map<string,vector<hit> >) );
      instance.SetNew(&new_maplEstringcOvectorlEhitgRsPgR);
      instance.SetNewArray(&newArray_maplEstringcOvectorlEhitgRsPgR);
      instance.SetDelete(&delete_maplEstringcOvectorlEhitgRsPgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOvectorlEhitgRsPgR);
      instance.SetDestructor(&destruct_maplEstringcOvectorlEhitgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,vector<hit> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<string,vector<hit> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOvectorlEhitgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<string,vector<hit> >*)0x0)->GetClass();
      maplEstringcOvectorlEhitgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOvectorlEhitgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOvectorlEhitgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,vector<hit> > : new map<string,vector<hit> >;
   }
   static void *newArray_maplEstringcOvectorlEhitgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,vector<hit> >[nElements] : new map<string,vector<hit> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOvectorlEhitgRsPgR(void *p) {
      delete ((map<string,vector<hit> >*)p);
   }
   static void deleteArray_maplEstringcOvectorlEhitgRsPgR(void *p) {
      delete [] ((map<string,vector<hit> >*)p);
   }
   static void destruct_maplEstringcOvectorlEhitgRsPgR(void *p) {
      typedef map<string,vector<hit> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<string,vector<hit> >

namespace ROOT {
   static TClass *maplEintcOvectorlEintgRsPgR_Dictionary();
   static void maplEintcOvectorlEintgRsPgR_TClassManip(TClass*);
   static void *new_maplEintcOvectorlEintgRsPgR(void *p = 0);
   static void *newArray_maplEintcOvectorlEintgRsPgR(Long_t size, void *p);
   static void delete_maplEintcOvectorlEintgRsPgR(void *p);
   static void deleteArray_maplEintcOvectorlEintgRsPgR(void *p);
   static void destruct_maplEintcOvectorlEintgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,vector<int> >*)
   {
      map<int,vector<int> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,vector<int> >));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,vector<int> >", -2, "map", 96,
                  typeid(map<int,vector<int> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOvectorlEintgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(map<int,vector<int> >) );
      instance.SetNew(&new_maplEintcOvectorlEintgRsPgR);
      instance.SetNewArray(&newArray_maplEintcOvectorlEintgRsPgR);
      instance.SetDelete(&delete_maplEintcOvectorlEintgRsPgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOvectorlEintgRsPgR);
      instance.SetDestructor(&destruct_maplEintcOvectorlEintgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,vector<int> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,vector<int> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOvectorlEintgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,vector<int> >*)0x0)->GetClass();
      maplEintcOvectorlEintgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOvectorlEintgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOvectorlEintgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,vector<int> > : new map<int,vector<int> >;
   }
   static void *newArray_maplEintcOvectorlEintgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,vector<int> >[nElements] : new map<int,vector<int> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOvectorlEintgRsPgR(void *p) {
      delete ((map<int,vector<int> >*)p);
   }
   static void deleteArray_maplEintcOvectorlEintgRsPgR(void *p) {
      delete [] ((map<int,vector<int> >*)p);
   }
   static void destruct_maplEintcOvectorlEintgRsPgR(void *p) {
      typedef map<int,vector<int> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,vector<int> >

namespace ROOT {
   static TClass *maplEintcOvectorlEdoublegRsPgR_Dictionary();
   static void maplEintcOvectorlEdoublegRsPgR_TClassManip(TClass*);
   static void *new_maplEintcOvectorlEdoublegRsPgR(void *p = 0);
   static void *newArray_maplEintcOvectorlEdoublegRsPgR(Long_t size, void *p);
   static void delete_maplEintcOvectorlEdoublegRsPgR(void *p);
   static void deleteArray_maplEintcOvectorlEdoublegRsPgR(void *p);
   static void destruct_maplEintcOvectorlEdoublegRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,vector<double> >*)
   {
      map<int,vector<double> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,vector<double> >));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,vector<double> >", -2, "map", 96,
                  typeid(map<int,vector<double> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOvectorlEdoublegRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(map<int,vector<double> >) );
      instance.SetNew(&new_maplEintcOvectorlEdoublegRsPgR);
      instance.SetNewArray(&newArray_maplEintcOvectorlEdoublegRsPgR);
      instance.SetDelete(&delete_maplEintcOvectorlEdoublegRsPgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOvectorlEdoublegRsPgR);
      instance.SetDestructor(&destruct_maplEintcOvectorlEdoublegRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,vector<double> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,vector<double> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOvectorlEdoublegRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,vector<double> >*)0x0)->GetClass();
      maplEintcOvectorlEdoublegRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOvectorlEdoublegRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOvectorlEdoublegRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,vector<double> > : new map<int,vector<double> >;
   }
   static void *newArray_maplEintcOvectorlEdoublegRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,vector<double> >[nElements] : new map<int,vector<double> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOvectorlEdoublegRsPgR(void *p) {
      delete ((map<int,vector<double> >*)p);
   }
   static void deleteArray_maplEintcOvectorlEdoublegRsPgR(void *p) {
      delete [] ((map<int,vector<double> >*)p);
   }
   static void destruct_maplEintcOvectorlEdoublegRsPgR(void *p) {
      typedef map<int,vector<double> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,vector<double> >

namespace ROOT {
   static TClass *maplEintcOdoublegR_Dictionary();
   static void maplEintcOdoublegR_TClassManip(TClass*);
   static void *new_maplEintcOdoublegR(void *p = 0);
   static void *newArray_maplEintcOdoublegR(Long_t size, void *p);
   static void delete_maplEintcOdoublegR(void *p);
   static void deleteArray_maplEintcOdoublegR(void *p);
   static void destruct_maplEintcOdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,double>*)
   {
      map<int,double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,double>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,double>", -2, "map", 96,
                  typeid(map<int,double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOdoublegR_Dictionary, isa_proxy, 4,
                  sizeof(map<int,double>) );
      instance.SetNew(&new_maplEintcOdoublegR);
      instance.SetNewArray(&newArray_maplEintcOdoublegR);
      instance.SetDelete(&delete_maplEintcOdoublegR);
      instance.SetDeleteArray(&deleteArray_maplEintcOdoublegR);
      instance.SetDestructor(&destruct_maplEintcOdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,double>*)0x0)->GetClass();
      maplEintcOdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,double> : new map<int,double>;
   }
   static void *newArray_maplEintcOdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,double>[nElements] : new map<int,double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOdoublegR(void *p) {
      delete ((map<int,double>*)p);
   }
   static void deleteArray_maplEintcOdoublegR(void *p) {
      delete [] ((map<int,double>*)p);
   }
   static void destruct_maplEintcOdoublegR(void *p) {
      typedef map<int,double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,double>

namespace {
  void TriggerDictionaryInitialization_struct_Impl() {
    static const char* headers[] = {
"struct.h",
0
    };
    static const char* includePaths[] = {
"/wd/sw/ROOT/root-6.16.00/root-6.16.00.binary/include",
"/wd/dune-it/enurec/analysis/kloe-simu/digitization/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "struct dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
struct __attribute__((annotate("$clingAutoload$struct.h")))  cell;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
struct __attribute__((annotate("$clingAutoload$struct.h")))  digit;
struct __attribute__((annotate("$clingAutoload$struct.h")))  particle;
struct __attribute__((annotate("$clingAutoload$struct.h")))  track;
struct __attribute__((annotate("$clingAutoload$struct.h")))  cluster;
namespace std{template <class _CharT> struct __attribute__((annotate("$clingAutoload$bits/char_traits.h")))  __attribute__((annotate("$clingAutoload$string")))  char_traits;
}
namespace std{template <typename _Tp> struct __attribute__((annotate("$clingAutoload$bits/stl_function.h")))  __attribute__((annotate("$clingAutoload$string")))  less;
}
namespace std{template <class _T1, class _T2> struct __attribute__((annotate("$clingAutoload$bits/stl_pair.h")))  __attribute__((annotate("$clingAutoload$string")))  pair;
}
struct __attribute__((annotate("$clingAutoload$struct.h")))  event;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "struct dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "struct.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"cell", payloadCode, "@",
"cluster", payloadCode, "@",
"digit", payloadCode, "@",
"event", payloadCode, "@",
"particle", payloadCode, "@",
"track", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("struct",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_struct_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_struct_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_struct() {
  TriggerDictionaryInitialization_struct_Impl();
}
