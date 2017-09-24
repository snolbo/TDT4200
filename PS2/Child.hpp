#pragma once

class Child : public Person {
private:

public:
  Child(std::string name, int age, int shoesize) : Person(name, age, shoesize) {
    // gjør noe
  }
  void play();
}
