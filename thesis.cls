\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{thesis}[2020/03/23, v1.0]

\RequirePackage[english]{babel}
\RequirePackage{calc}
\RequirePackage{etoolbox}
\RequirePackage{regexpatch}
\RequirePackage{amsmath}
\RequirePackage{unicode-math}
\RequirePackage{mathtools}
\RequirePackage{graphicx}
\RequirePackage{siunitx}
\RequirePackage{setspace}

\LoadClass[a4paper, 12pt, twoside]{memoir}

\setmainfont[Path = /usr/share/fonts/TTF/,
  Ligatures = TeX,
  Extension = .ttf,
  UprightFont = *-Light,
  BoldFont = *-Medium,
  ItalicFont = *-Italic,
  BoldItalicFont = *-MediumItalic]
  {IBMPlexSerif}
\setmonofont[Scale = MatchLowercase,
	Extension = .ttf]{RobotoMono-Light}

\RequirePackage[left=3.5cm,%
                marginparwidth=4cm,%
                marginparsep=0.5cm,%
                right=1.3cm,%
                top=2cm,%
                bottom=2.5cm,%
                headheight=21pt,%
                headsep=2\baselineskip,%
                columnsep=2em,%
                letterpaper]{geometry}%

\RequirePackage[labelfont={bf},%
                labelsep=period,%
                justification=raggedright,%
                singlelinecheck=false,%
                tableposition=top,font=small]
                {caption}

\bibliographystyle{vancouver}
