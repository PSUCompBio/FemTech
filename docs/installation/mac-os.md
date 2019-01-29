---
description: How to build FemTech on Mac OS
---

# Mac OS

 Download and install Xcode with developer command line tools \([https://www.embarcadero.com/starthere/berlin/mobdevsetup/ios/en/installing\_the\_xcode\_command\_line\_tools\_on\_a\_mac.html](https://www.embarcadero.com/starthere/berlin/mobdevsetup/ios/en/installing_the_xcode_command_line_tools_on_a_mac.html)\)

 CMake \([https://cmake.org/download/](https://cmake.org/download/)\)

 open terminal \(navigate to desired directory\)

```text
git clone https://github.com/PSUCompBio/FemTech
cd FemTech
mkdir build
cd bulid
ccmake .. (configure as needed)
make -j8
```



