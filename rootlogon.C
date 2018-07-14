void rootlogon()
{
  gInterpreter->LoadMacro("loader/loader.C+");
  gInterpreter->LoadMacro("CalDigi/CalDigi.cpp+");
  gInterpreter->LoadMacro("STTDigi/STTDigi.cpp+");
}