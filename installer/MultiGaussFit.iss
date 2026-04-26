[Setup]
; NOTE: The value of AppId uniquely identifies this application.
AppId={{A3D1E0C9-2B5C-4DF7-9E45-8C1A2F5EABC1}
AppName=MultiGaussFit
AppVersion=1.0
AppPublisher=JELC
DefaultDirName={autopf}\MultiGaussFit
DisableProgramGroupPage=yes
; Output folder for the installer
OutputDir=..\dist
OutputBaseFilename=MultiGaussFit_Setup
Compression=lzma
SolidCompression=yes
WizardStyle=modern
SetupIconFile=..\logo.ico

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"
Name: "spanish"; MessagesFile: "compiler:Languages\Spanish.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Files]
Source: "..\dist\MultiGaussFit\MultiGaussFit.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\dist\MultiGaussFit\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs
; NOTE: Don't use "Flags: ignoreversion" on any shared system files

[Icons]
Name: "{autoprograms}\MultiGaussFit"; Filename: "{app}\MultiGaussFit.exe"
Name: "{autodesktop}\MultiGaussFit"; Filename: "{app}\MultiGaussFit.exe"; Tasks: desktopicon

[Run]
Filename: "{app}\MultiGaussFit.exe"; Description: "{cm:LaunchProgram,MultiGaussFit}"; Flags: nowait postinstall skipifsilent
