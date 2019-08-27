//###################################################################
//### Name:        ShowerElementHolder                           ### 
//### Author:      Dominic Barker                                 ###
//### Date:        15.07.19                                       ###
//### Description: Class to holder the standard shower property   ###
//###              information. Used in SBNShower_modle and       ###
//###              corresponding tools.                           ###
//###################################################################

#ifndef ShowerElementHolder_HH
#define ShowerElementHolder_HH

//Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

//C++ Inlcudes 
#include <iostream> 
#include <map>
#include <string> 
#include <memory>
#include <iomanip>
#include <cxxabi.h>

namespace reco {
  namespace shower {
    class ShowerElementBase;
    template <class T> class ShowerElementAccessor;     
    template <class T> class ShowerDataProduct; 
    template <class T, class T2> class ShowerProperty;
    class ShowerElementHolder;
  }
}

class reco::shower::ShowerElementBase {


public:

  virtual ~ShowerElementBase() noexcept = default;
  
  virtual bool CheckTag(){
    throw cet::exception("ShowerElementHolder") << "Trying to check an element that is not a product" << std::endl;
  }
  virtual void SetCheckTag(bool& check){
    throw cet::exception("ShowerElementHolder") << "Trying to set an element that is not a product" << std::endl;
  }

  virtual std::string GetType() = 0;

  //Check if the element has been set.
  bool CheckShowerElement(){
    if(elementPtr) return true;
    else return false;
  }

   void Clear(){
    elementPtr    = 0;
  }


protected:

  bool elementPtr;
  
};

//This is a template class which holds a shower property. This holds any object e.g. std::vector<double>, double, TVector3 
//and holds various information which helps the showerproperty holder access the elements. A user should not require any part 
//of this class.  
template <class T>
class reco::shower::ShowerElementAccessor : public reco::shower::ShowerElementBase {

public: 

  //Set the element in the holder
  void SetShowerElement(T& Element){
    element = Element;
    this->elementPtr = 1;
  }
  
  //Fill Element with the element that the holder holds.
  int GetShowerElement(T& Element){
    if(this->elementPtr){
      Element = element;
      return 0;
    }
    else{
      return 1;
    }
  }

  //Return a copy of the shower element. 
  T GetShowerElement(){
    if(this->elementPtr){
      return element;
    }
    else{
      return T();
    }
  } 

  //Return the type as a string.
  std::string GetType() override {
    int status = -9;
    return abi::__cxa_demangle(typeid(element).name(),NULL,NULL,&status);  
  }

protected:
  T   element; 
};

//This class holds shower data products which have the potential to be saved in the art::Event e.g. recob::Track. Note the product itself must be store in the element holder as the object will be destoryed in the CalculateProperty Section otherwise. Associtations can be made during Calculate property tool stage. 
template <class T> 
class reco::shower::ShowerDataProduct : public reco::shower::ShowerElementAccessor<T>{

public:
  
  ShowerDataProduct(){
    this->elementPtr      = 0;
    checktag            = false;
  }

  ShowerDataProduct(T& Element, bool& Checktag){
    this->elementPtr      = 1;
    this->element         = Element;
    checktag            = Checktag;
  }


  void Clear(){
    this->element = T();
    this->elementPtr    = 0;
  }

  //Check if we should check the dataproduct in the end.
  bool CheckTag(){
    return checktag;
  }

  //Set if we should check the data product in the end.
  void SetCheckTag(bool& Checktag){
    checktag = Checktag;
  }

  private:
  bool checktag;
};

//This class holds shower properties e.g. ShowerDirection. The user must define the associated error  
template <class T, class T2>
class reco::shower::ShowerProperty : public reco::shower::ShowerElementAccessor<T>{

public:

  ShowerProperty(){
    this->elementPtr = 0;
  }
  
  ShowerProperty(T& Element, T2& ElementErr){
    this->elementPtr = 1;
    this->element    = Element;
    propertyErr      = ElementErr;
  }
  
  //Fill the property error as long as it has been set.
  int GetShowerPropertyError(T2& ElementErr){
    if(this->elementPtr){
      ElementErr = propertyErr;
      return 0;
    }
    else{
      return 1;
    }
  }

  //Set the properties. Note you cannot set an property without an error.
  void SetShowerProperty(T& Element, T2& ElementErr){
    this->element    = Element;
    this->elementPtr = 1;
    propertyErr      = ElementErr;
  }

  void Clear(){
    this->element = T();
    this->elementPtr = 0;
  }
  
private:
  T2   propertyErr;

};


//Class to holder all the reco::shower::ShowerElement objects. This is essentially a map from a string the object so people can 
//add an object in a tool and get it back later. 
class reco::shower::ShowerElementHolder{

public:
  
  //Getter function for accessing the shower property e..g the direction ShowerElementHolder.GetElement("MyShowerValue"); The name is used access the value and precise names are required for a complete shower in sbnshower: ShowerStartPosition, ShowerDirection, ShowerEnergy ,ShowerdEdx.
  template <class T >
  int GetElement(std::string Name, T& Element){
    if(showerproperties.find(Name) != showerproperties.end()){
      if(showerproperties[Name]->CheckShowerElement()){
	  reco::shower::ShowerElementAccessor<T> *showerprop = dynamic_cast<reco::shower::ShowerElementAccessor<T> *>(showerproperties[Name].get());
	  if(showerprop == NULL){
	    throw cet::exception("ShowerElementHolder") << "Trying to get Element: " << Name << ". This element you are filling is not the correct type" << std::endl;
	  }
	  showerprop->GetShowerElement(Element);
	  return 0;
	}
      else{
	mf::LogWarning("ShowerElementHolder") << "Trying to get Element " << Name << ". This elment has not been filled" << std::endl;
	return 1;
      }
    }
    else if(showerdataproducts.find(Name) != showerdataproducts.end()){
      if(showerdataproducts[Name]->CheckShowerElement()){
	reco::shower::ShowerElementAccessor<T> *showerprop = dynamic_cast<reco::shower::ShowerElementAccessor<T> *>(showerdataproducts[Name].get());
	if(showerprop == NULL){
	  throw cet::exception("ShowerElementHolder") << "Trying to get Element: " << Name << ". This element you are filling is not the correct type" << std::endl;
	}
	showerprop->GetShowerElement(Element);
	return 0;
      }
      else{
	mf::LogWarning("ShowerElementHolder") << "Trying to get Element " << Name << ". This elment has not been filled" << std::endl;
	return 1;
      }
    }
    throw cet::exception("ShowerElementHolder") << "Trying to get Element: " << Name << ". This element does not exist in the element holder" << std::endl;
    return 1;
  }

  //Alternative get function that returns the object. Not recommended.
  template <class T > 
  T GetElement(std::string Name){
    if(showerproperties.find(Name) != showerproperties.end()){ 
      if(showerproperties[Name]->CheckShowerElement()){
	reco::shower::ShowerElementAccessor<T> *showerprop = dynamic_cast<reco::shower::ShowerElementAccessor<T> *>(showerproperties[Name].get());
	if(showerprop == NULL){
	  throw cet::exception("ShowerElementHolder") << "Trying to get Element: " << Name << ". This element you are filling is not the correct type" << std::endl;
	}
	return showerprop->GetShowerElement();
      }
    }
    else if(showerdataproducts.find(Name) != showerdataproducts.end()){   
      if(showerdataproducts[Name]->CheckShowerElement()){ 
	reco::shower::ShowerElementAccessor<T> *showerprop = dynamic_cast<reco::shower::ShowerElementAccessor<T> *>(showerdataproducts[Name].get());
	if(showerprop == NULL){
	  throw cet::exception("ShowerElementHolder") << "Trying to get Element: " << Name << ". This element you are filling is not the correct type" << std::endl;
	}
	return showerprop->GetShowerElement();
      }
    }	
    throw cet::exception("ShowerElementHolder") << "Trying to get Element: " << Name << ". This element does not exist in the element holder" << std::endl;
    return 1;
  }
  
  //Getter function for accessing the shower property error e.g the direction ShowerElementHolder.GetElement("MyShowerValue");
  template <class T, class T2>
  int GetElementAndError(std::string Name, T& Element,  T2& ElementErr){
    if(showerproperties.find(Name) == showerproperties.end()){
      mf::LogError("ShowerElementHolder") << "Trying to get Element Error: " << Name << ". This elment does not exist in the element holder" << std::endl;
      return 1;
    }
    reco::shower::ShowerProperty<T,T2> *showerprop = dynamic_cast<reco::shower::ShowerProperty<T,T2> *>(showerproperties[Name].get());
    showerprop->GetShowerElement(Element);  
    showerprop->GetShowerPropertyError(ElementErr);
    return 0;
  }


  //This sets the value of the data product. Just give a name and a object
  //e.g. TVector3 ShowerElementHolder.SetElement((TVector3) StartPosition, "StartPosition");
  template <class T>
  void SetElement(T& dataproduct, std::string Name, bool checktag=false){
    
    if(showerdataproducts.find(Name) != showerdataproducts.end()){
      reco::shower::ShowerDataProduct<T>* showerdataprod = dynamic_cast<reco::shower::ShowerDataProduct<T> *>(showerdataproducts[Name].get());
      showerdataprod->SetShowerElement(dataproduct);
      showerdataprod->SetCheckTag(checktag);
      return;
    }
    else{
      showerdataproducts[Name] = std::unique_ptr<reco::shower::ShowerDataProduct<T> >(new reco::shower::ShowerDataProduct<T>(dataproduct,checktag));
      return;
    }
  }

  //This sets the value of the property. Just give a name and a object
  //e.g. TVector3 ShowerElementHolder.SetElement((art::Ptr<recob::Track>) track, "StartPosition", save);
  template <class T, class T2>
  void SetElement(T& propertyval, T2& propertyvalerror, std::string Name){

    if(showerproperties.find(Name) != showerproperties.end()){
      reco::shower::ShowerProperty<T,T2>* showerprop = dynamic_cast<reco::shower::ShowerProperty<T,T2> *>(showerproperties[Name].get());
      showerprop->SetShowerProperty(propertyval,propertyvalerror);
    return;
    }
    else{
      showerproperties[Name] = std::unique_ptr<reco::shower::ShowerProperty<T,T2> >(new reco::shower::ShowerProperty<T,T2>(propertyval,propertyvalerror));
      return;
    }
  }

  //Check that a property is filled 
  bool CheckElement(std::string Name){
    if(showerproperties.find(Name) != showerproperties.end()){
      return showerproperties[Name]->CheckShowerElement();
    }
    if(showerdataproducts.find(Name) != showerdataproducts.end()){
      return showerdataproducts[Name]->CheckShowerElement();
    }
    return false;
  }

  //Check All the properties
  bool CheckAllElements(){
    bool checked = true;
    for(auto const& showerprop: showerproperties){
     checked *= showerprop.second->CheckShowerElement();
    }
    for(auto const& showerdataprod: showerdataproducts){
      checked *= showerdataprod.second->CheckShowerElement();
    }
    return checked;
  }

  
  //Clear Fucntion. This does not delete the element.
  void ClearElement(std::string Name){
    if(showerproperties.find(Name) != showerproperties.end()){
      return showerproperties[Name]->Clear();
    }
    if(showerdataproducts.find(Name) != showerdataproducts.end()){
      return showerdataproducts[Name]->Clear();
    }
    mf::LogError("ShowerElementHolder") << "Trying to clear Element: " << Name << ". This element does not exist in the element holder" << std::endl;
    return;
  }

  //Clear all the shower properties. This does not delete the element.
  void ClearAll(){
    for(auto const& showerprop: showerproperties){
      (showerprop.second)->Clear();
    }
    for(auto const& showerdataproducts: showerdataproducts){
      (showerdataproducts.second)->Clear();
    }
  }

  //Find if the product is one what is being stored.
  bool CheckElementTag(std::string Name){
    if(showerdataproducts.find(Name) != showerdataproducts.end()){
    return showerdataproducts[Name]->CheckTag();
    }
    return false;
  }
 
  //Delete a product. I see no reason for it.
  void DeleteElement(std::string Name){
    if(showerdataproducts.find(Name) != showerdataproducts.end()){
      showerdataproducts[Name].reset(nullptr);
      return;
    }
    if(showerproperties.find(Name) != showerproperties.end()){
      showerproperties[Name].reset(nullptr);
      return;
    }
    mf::LogError("ShowerElementHolder") << "Trying to delete Element: " << Name << ". This element does not exist in the element holder" << std::endl;
    return;
  }

  //Set the indicator saying if the shower is going to be stored.
  void SetElementTag(std::string Name, bool checkelement){
    if(showerdataproducts.find(Name) != showerdataproducts.end()){
      showerdataproducts[Name]->SetCheckTag(checkelement);
      return;
    }
    mf::LogError("ShowerElementHolder") << "Trying set the checking of the data product: " << Name << ". This data product does not exist in the element holder" << std::endl;
    return;
  }

  bool CheckAllElementTags(){
    bool checked = true;
    for(auto const& showerdataproduct: showerdataproducts){
      bool check  = showerdataproduct.second->CheckTag();
      if(check){
	bool elementset = showerdataproduct.second->CheckShowerElement();
	if(!elementset){
	  mf::LogError("ShowerElementHolder") << "The following element is not set and was asked to be checked: " << showerdataproduct.first << std::endl;
	  checked = false;
	}
      }
    }
    return checked;
  }

  //Set the shower number. This is required the association making.
  void SetShowerNumber(int& shower_iter){
    showernumber = shower_iter;
  }
  
  //Get the shower number.
  int GetShowerNumber(){
    return showernumber;
  }

  void PrintElement(std::string Name){
    if(showerdataproducts.find(Name) != showerdataproducts.end()){
      std::string Type = showerdataproducts[Name]->GetType();
      std::cout << "Element Name: " << Name << " Type: " << Type << std::endl;
      return;
    }
    if(showerproperties.find(Name) != showerproperties.end()){
      std::string Type = showerproperties[Name]->GetType();
      std::cout << "Element Name: " << Name << " Type: " << Type << std::endl;
      return;
    }
    mf::LogError("ShowerElementHolder") << "Trying to print Element: " << Name << ". This element does not exist in the element holder" << std::endl;
    return;
  }

  //This function will print out all the elements and there types for the user to check. 
  void PrintElements(){

    unsigned int maxname = 0;
    for(auto const& showerprop: showerproperties){
      if(showerprop.first.size() > maxname){
	maxname = showerprop.first.size();
      }
    }
    for(auto const& showerdataprod: showerdataproducts){
      if(showerdataprod.first.size() > maxname){
	maxname = showerdataprod.first.size();
      }
    }

    std::map<std::string,std::string> Type_showerprops;
    std::map<std::string,std::string> Type_showerdataprods;
    for(auto const& showerprop: showerproperties){
      std::string Type = (showerprop.second)->GetType();
      Type_showerprops[showerprop.first] = Type;
    }
    for(auto const& showerdataprod: showerdataproducts){
      std::string Type = (showerdataprod.second)->GetType();
      Type_showerdataprods[showerdataprod.first] = Type;
    }

    unsigned int maxtype = 0;
    for(auto const& Type_showerprop: Type_showerprops){
      if(Type_showerprop.second.size() > maxtype){
	maxtype = Type_showerprop.second.size();
      }
    }
    for(auto const& Type_showerdataprod: Type_showerdataprods){
      if(Type_showerdataprod.second.size() > maxtype){
	maxtype = Type_showerdataprod.second.size();
      }
    }

    unsigned int n = maxname + maxtype + 33;
    std::cout << std::left << std::setfill('*') << std::setw(n-1) << "*" <<std::endl;
    std::cout << "Elements in the element holder" << std::endl;
    std::cout << std::left << std::setfill('*') << std::setw(n-1) << "*" <<std::endl;
    for(auto const& Type_showerprop: Type_showerprops){
      std::cout << std::left << std::setfill(' ') << std::setw(21) << "* Property Name: " << std::setw(maxname) << Type_showerprop.first;
      std::cout << std::left << std::setfill(' ') << " * Type: " << std::setw(maxtype) << Type_showerprop.second <<  " * " << std::endl;
    }
    for(auto const& Type_showerdataprod: Type_showerdataprods){
      std::cout << std::left << std::setfill(' ') << std::setw(maxname) << std::setw(21)  << "* Data Product Name: " << std::setw(maxname) << Type_showerdataprod.first;
      std::cout << std::left << std::setfill(' ') << " * Type: " << std::setw(maxtype) <<  Type_showerdataprod.second << " *" << std::endl;
    }
    std::cout << std::left << std::setfill('*') << std::setw(n-1) << "*" <<std::endl;
    std::cout << std::setfill(' ');
    std::cout << std::setw(0); 
    return;
  }

private:

  //Storage for all the shower properties.
  std::map<std::string,std::unique_ptr<reco::shower::ShowerElementBase> > showerproperties;
  
  //Storage for all the data products
  std::map<std::string,std::unique_ptr<reco::shower::ShowerElementBase> > showerdataproducts;

  //Shower ID number. Use this to set ptr makers.
  int showernumber;

};


#endif
