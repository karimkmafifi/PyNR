
//============================================================================
// stl_ptr_GN.hpp -*- C++ -*-; pointer class
//
// Copyright (C) 2006, 2007, 2008, 2009, 2010 Gerd Neudert
//
// This file is part of fconv.
// fconv is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------------
//
// author:      Gerd Neudert
//              gneudert(place_at_here)web.de
// supervisor:  Prof. Dr. Gerhard Klebe
//              klebe(place_at_here)staff.uni-marburg.de
// Department of Pharmaceutical Chemistry, Philipps-University of Marburg
//
// This library is part of the program fconv. Please see the
// documentation in the file fconv.cpp first!
// This library implements a pointer object. This is NO smart pointer!
// Please do not implement operator T*(), as deletes for stl_ptr would then
// result in deletes of the objects they are pointing to (use method kill()
// to delete the pointed objects).
// The comparison operators call the comparison operators of the pointed
// objects and do not compare the pointers!
//============================================================================


#ifndef __STLPTRGN
#define __STLPTRGN

#include<sstream>

using namespace std;

//==================================================================================================
//!Forward-Deklarationen:
//==================================================================================================

template<class T> class stl_ptr;
template<class T> ostream &operator<<(ostream &os,stl_ptr<T> const& refptr);

//==================================================================================================
//!Deklaration der Klasse stl_ptr:
//==================================================================================================

template<class T>
class stl_ptr {
    private:
        T *ptr;
    public:
        stl_ptr(T *ptr = NULL);
        stl_ptr(T &ref);
        stl_ptr(stl_ptr<T> const& refptr);
        ~stl_ptr();
        
        inline void kill(); //loescht das Objekt, auf das ptr zeigt (delete loescht nur den stl_ptr)
        inline void zero_kill(); //Objekt loeschen und Zeiger auf Null setzen
        
        inline T* get_pointer() const;
        
        inline const stl_ptr<T>& operator=(stl_ptr<T> const& rechts); //Zuweisungsoperator
        inline const stl_ptr<T>& operator=(T &rechts); //Zuweisung ber das Objekt, auf das gezeigt werden soll
        inline const stl_ptr<T>& operator=(T *rechts); //direkte Zuweisung eines Zeigers
        T& operator*(); //Dereferenzierungsoperator (liefert eine Referenz auf das Objekt!!!)
        //!wird auch als const gebraucht (wenn 'this' const sein muss, z.B. fuer stl-algorithmen)
        T& operator*() const;
        T* operator->();
        T* operator->() const;
        inline bool equal_addr(stl_ptr<T> const& rechts);
        inline bool equal_addr(stl_ptr<T> const& rechts) const;
        inline bool operator==(stl_ptr<T> const& rechts) const; //Vergleicht die Objekte auf die gezeigt wird
        inline bool operator<(stl_ptr<T> const& rechts) const; //Die Vergleichsoperatoren muessen also fuer T definiert sein!!!
        inline bool operator>(stl_ptr<T> const& rechts) const;
        inline bool operator!=(stl_ptr<T> const& rechts) const;
        inline bool zero() const; //prueft auf Nullzeiger
        friend ostream &operator<< <>(ostream &os,stl_ptr<T> const& refptr);
};


//==================================================================================================
//!Definitionen der Klasse stl_ptr:
//==================================================================================================

template<class T> 
stl_ptr<T>::stl_ptr(T *p):ptr(p) {}

template<class T> 
stl_ptr<T>::stl_ptr(T &ref):ptr(&ref) {}

template<class T> 
stl_ptr<T>::stl_ptr(const stl_ptr<T> &refptr):ptr(refptr.ptr) {} //!echte Kopie des Zeigers => 2 Zeiger auf das gleiche Objekt!!!

template<class T> 
stl_ptr<T>::~stl_ptr() {} //!Destruktor ohne automatisches delete (Objekt hinter ptr wird nicht zerstoert!!!)

template<class T> 
void stl_ptr<T>::kill() {delete ptr;}

template<class T> 
void stl_ptr<T>::zero_kill() {delete ptr; ptr = NULL;}

template<class T> 
T* stl_ptr<T>::get_pointer() const {return ptr;}

template<class T> 
const stl_ptr<T>& stl_ptr<T>::operator=(stl_ptr<T> const& rechts) {
    if (this == &rechts) return *this;
    ptr = rechts.ptr; //!echte Kopie des Zeigers => 2 Zeiger auf das gleiche Objekt
    return *this;
}

template<class T> 
const stl_ptr<T>& stl_ptr<T>::operator=(T &rechts) {
    ptr = &rechts;
    return *this;
}

template<class T> 
const stl_ptr<T>& stl_ptr<T>::operator=(T *rechts) {
    ptr = rechts;
    return *this;
}

template<class T> 
T& stl_ptr<T>::operator*() {
    return *ptr;
}

template<class T> 
T& stl_ptr<T>::operator*() const {
    return *ptr;
}

template<class T> 
T* stl_ptr<T>::operator->() {
    return ptr;
}

template<class T> 
T* stl_ptr<T>::operator->() const {
    return ptr;
}

template<class T>
bool stl_ptr<T>::equal_addr(stl_ptr<T> const& rechts) {
    if (ptr == rechts.ptr) return true;
    return false;
}

template<class T>
bool stl_ptr<T>::equal_addr(stl_ptr<T> const& rechts) const {
    if (ptr == rechts.ptr) return true;
    return false;
}

template<class T> 
bool stl_ptr<T>::operator==(stl_ptr<T> const& rechts) const {
    if (*ptr == *rechts) return true;
    return false;
}

template<class T> 
bool stl_ptr<T>::operator!=(stl_ptr<T> const& rechts) const {
    if (*ptr != *rechts) return true;
    return false;
}

template<class T> 
bool stl_ptr<T>::operator<(stl_ptr<T> const& rechts) const {
    if (*ptr < *rechts) return true;
    return false;
}

template<class T> 
bool stl_ptr<T>::operator>(stl_ptr<T> const& rechts) const {
    if (*ptr > *rechts) return true;
    return false;
}

template<class T> 
bool stl_ptr<T>::zero() const {
    if (ptr == NULL) return true;
    return false;
}

template<class T> 
ostream &operator<<(ostream &os,stl_ptr<T> const& refptr) {
    os << (*refptr.ptr);
    return os;
}

#endif //__STLPTRGN
