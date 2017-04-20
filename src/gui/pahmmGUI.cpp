// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------


// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
    #pragma hdrstop
#endif

// for all others, include the necessary headers (this file is usually all you
// need because it includes almost all "standard" wxWidgets headers)
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include "wx/filename.h"

#if wxUSE_FILEDLG
    #include "wx/filedlg.h"
#endif // wxUSE_FILEDLG

#if USE_FILEDLG_GENERIC
    #include "wx/generic/filedlgg.h"
#endif // USE_FILEDLG_GENERIC

#include "wx/dataview.h"
#include "wx/datetime.h"
#include "wx/splitter.h"
#include "wx/aboutdlg.h"
#include "wx/colordlg.h"
#include "wx/choicdlg.h"
#include "wx/numdlg.h"
#include "wx/spinctrl.h"
#include "wx/imaglist.h"
#include "wx/itemattr.h"
#include "wx/notebook.h"
#include "wx/hashmap.h"


#include "gui/pahmmGUI.h"

// ----------------------------------------------------------------------------
// event tables and other macros for wxWidgets
// ----------------------------------------------------------------------------

// the event tables connect the wxWidgets events with the functions (event
// handlers) which process them. It can be also done at run-time, but for the
// simple menu events like this the static method is much simpler.
wxBEGIN_EVENT_TABLE(PaHmmFrame, wxFrame)
    EVT_MENU(PaHmm_Quit,  PaHmmFrame::OnQuit)
    EVT_MENU(PaHmm_About, PaHmmFrame::OnAbout)
    
    EVT_MENU(PaHmm_FileOpenSingle, PaHmmFrame::FileOpen)
    EVT_MENU(PaHmm_FileOpenMultiple, PaHmmFrame::FilesOpen)

    EVT_MENU(PaHmm_GTRModelSelect, PaHmmFrame::OnGTRModelSelect)
    EVT_MENU(PaHmm_HKYModelSelect, PaHmmFrame::OnHKYModelSelect)
    EVT_MENU(PaHmm_LGModelSelect, PaHmmFrame::OnLGModelSelect)
    EVT_MENU(PaHmm_GammaSelect, PaHmmFrame::OnGammaSelect)
    EVT_MENU(PaHmm_ModelParamsSelect, PaHmmFrame::OnModelParams)
    EVT_MENU(PaHmm_IndelParamsSelect, PaHmmFrame::OnIndelParams)
wxEND_EVENT_TABLE()

// Create a new application object: this macro will allow wxWidgets to create
// the application object during program execution (it's better than using a
// static object for many reasons) and also implements the accessor function
// wxGetApp() which will return the reference of the right type (i.e. PaHmmApp and
// not wxApp)
wxIMPLEMENT_APP(PaHmmApp);


// 'Main program' equivalent: the program execution "starts" here
bool PaHmmApp::OnInit()
{
    // call the base class initialization method, currently it only parses a
    // few common command-line options but it could be do more in the future
    if ( !wxApp::OnInit() )
        return false;

    // create the main application window
    PaHmmFrame *frame = new PaHmmFrame("paHMM-Tree-GUI");

    // and show it (the frames, unlike simple controls, are not shown when
    // created initially)
    frame->Show(true);

    // success: wxApp::OnRun() will be called which will enter the main message
    // loop and the application will run. If we returned false here, the
    // application would exit immediately.
    return true;
}

// ----------------------------------------------------------------------------
// main frame
// ----------------------------------------------------------------------------

// frame constructor
PaHmmFrame::PaHmmFrame(const wxString& title)
       : wxFrame(NULL, wxID_ANY, title)
{
    // create a menu bar
    wxMenu *fileMenu = new wxMenu;
    wxMenu *modelMenu = new wxMenu;

    // the "About" item should be in the help menu
    wxMenu *helpMenu = new wxMenu;
    helpMenu->Append(PaHmm_About, "&About\tF1", "Show about dialog");


    wxTextCtrl    *m_log;
    wxLog         *m_logOld;

    m_log = new wxTextCtrl( this, wxID_ANY, wxT("This is the log window.\n"),
                            wxPoint(5,260), wxSize(630,100),
                            wxTE_MULTILINE | wxTE_READONLY);

    m_logOld = wxLog::SetActiveTarget( new wxLogTextCtrl( m_log ) );

    fileMenu->Append(PaHmm_FileOpenSingle,  wxT("&Open file\tCtrl-O"));
    fileMenu->Append(PaHmm_FileOpenMultiple,  wxT("Open &files\tShift-Ctrl-O"));
    fileMenu->Append(PaHmm_Quit, "E&xit\tAlt-X", "Quit this program");

    modelMenu->AppendRadioItem(PaHmm_GTRModelSelect, wxT("GTR Model (nucleotides)"));
    modelMenu->AppendRadioItem(PaHmm_HKYModelSelect, wxT("HKY85 Model (nucleotides)"));
    modelMenu->AppendRadioItem(PaHmm_LGModelSelect, wxT("LG Model (amino acids)"));
    modelMenu->AppendSeparator();
    modelMenu->AppendCheckItem(PaHmm_GammaSelect, wxT("Gamma rate heterogeneity"));
    modelMenu->AppendSeparator();
    modelMenu->Append(PaHmm_ModelParamsSelect, wxT("Substitution Model Parameters"));
    modelMenu->Append(PaHmm_IndelParamsSelect, wxT("Indel Model Parameters"));

    // now append the freshly created menu to the menu bar...
    wxMenuBar *menuBar = new wxMenuBar();
    menuBar->Append(fileMenu, "&File");
    menuBar->Append(modelMenu, "&Model");
    menuBar->Append(helpMenu, "&Help");

    // ... and attach this menu bar to the frame
    SetMenuBar(menuBar);
    // If menus are not available add a button to access the about box
    wxSizer* sizer = new wxBoxSizer(wxVERTICAL);
    //wxButton* aboutBtn = new wxButton(this, wxID_ANY, "About...");
    //aboutBtn->Bind(wxEVT_BUTTON, &PaHmmFrame::OnAbout, this);
    //sizer->Add(aboutBtn, wxSizerFlags().Center());
    wxDataViewListCtrl* lc =
                new wxDataViewListCtrl( this, wxID_ANY, wxDefaultPosition,
                                        wxDefaultSize);

            MyListStoreDerivedModel* page2_model = new MyListStoreDerivedModel();
            lc->AssociateModel(page2_model);
            page2_model->DecRef();

            lc->AppendToggleColumn( "Toggle" );
            lc->AppendTextColumn( "Text" );
            lc->AppendProgressColumn( "Progress" );

            wxVector<wxVariant> data;
            for (unsigned int i=0; i<10; i++)
            {
                data.clear();
                data.push_back( (i%3) == 0 );
                data.push_back( wxString::Format("row %d", i) );
                data.push_back( long(5*i) );

                lc->AppendItem( data );
            }

    sizer->Add( lc, 2, wxEXPAND | wxALL , 5 );
    sizer->Add( m_log, 1, wxALL | wxEXPAND, 5 );

    SetSizerAndFit(sizer);
}

void PaHmmFrame::FileOpen(wxCommandEvent& WXUNUSED(event) )
{
    wxFileDialog dialog
                 (
                    this,
                    wxT("Select FASTA file to open"),
                    wxEmptyString,
                    wxEmptyString,
                    wxT("FASTA files (*.fas)|*.fasta")
                 );

    dialog.CentreOnParent();
    dialog.SetDirectory(wxGetHomeDir());

    if (dialog.ShowModal() == wxID_OK)
    {
        wxString info;
        wxWindow * const extra = dialog.GetExtraControl();
        info.Printf(wxT("Full file name: %s\n")
                    wxT("Path: %s\n")
                    wxT("Name: %s\n"),
                    dialog.GetPath().c_str(),
                    dialog.GetDirectory().c_str(),
                    dialog.GetFilename().c_str());
        wxMessageDialog dialog2(this, info, wxT("Selected file"));
        dialog2.ShowModal();
    }
}


void PaHmmFrame::FilesOpen(wxCommandEvent& WXUNUSED(event) )
{
    wxString wildcards = wxT("FASTA files (*.fas)|*.fasta");
    wxFileDialog dialog(this, wxT("Open multiple FASTA files"),
                        wxEmptyString, wxEmptyString, wildcards,
                        wxFD_OPEN|wxFD_MULTIPLE);

    if (dialog.ShowModal() == wxID_OK)
    {
        wxArrayString paths, filenames;

        dialog.GetPaths(paths);
        dialog.GetFilenames(filenames);

        wxString msg, s;
        size_t count = paths.GetCount();
        for ( size_t n = 0; n < count; n++ )
        {
            s.Printf(wxT("File %d: %s (%s)\n"),
                     (int)n, paths[n].c_str(), filenames[n].c_str());

            msg += s;
        }
        s.Printf(wxT("Filter index: %d"), dialog.GetFilterIndex());
        msg += s;

        wxMessageDialog dialog2(this, msg, wxT("Selected files"));
        dialog2.ShowModal();
    }
}



// event handlers

void PaHmmFrame::OnGTRModelSelect(wxCommandEvent& WXUNUSED(event))
{
}
void PaHmmFrame::OnHKYModelSelect(wxCommandEvent& WXUNUSED(event))
{
}
void PaHmmFrame::OnLGModelSelect(wxCommandEvent& WXUNUSED(event))
{
}
void PaHmmFrame::OnGammaSelect(wxCommandEvent& WXUNUSED(event))
{
}

void PaHmmFrame::OnModelParams(wxCommandEvent& WXUNUSED(event))
{
}
void PaHmmFrame::OnIndelParams(wxCommandEvent& WXUNUSED(event))
{
}

void PaHmmFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
{
    // true is to force the frame to close
    Close(true);
}

void PaHmmFrame::OnAbout(wxCommandEvent& WXUNUSED(event))
{
    wxMessageBox(wxString::Format
                 (
                    "This is the GUI for paHMM-Tree - a tree inferene program\n"
                    "running under %s.",
                    wxGetOsDescription()
                 ),
                 "About paHMM-Tree-GUI",
                 wxOK | wxICON_INFORMATION,
                 this);
}


bool MyListStoreDerivedModel::IsEnabledByRow(unsigned int row, unsigned int col) const
{
    // disabled the last two checkboxes
    return !(col == 0 && 8 <= row && row <= 9);
}
