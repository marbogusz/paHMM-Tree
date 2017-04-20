#ifndef __GUI__
#define __GUI__

#ifdef __WXUNIVERSAL__
    #define USE_WXUNIVERSAL 1
#else
    #define USE_WXUNIVERSAL 0
#endif

#ifdef WXUSINGDLL
    #define USE_DLL 1
#else
    #define USE_DLL 0
#endif

#if defined(__WXMSW__)
    #define USE_WXMSW 1
#else
    #define USE_WXMSW 0
#endif

#ifdef __WXMAC__
    #define USE_WXMAC 1
#else
    #define USE_WXMAC 0
#endif

#if USE_NATIVE_FONT_DIALOG_FOR_MACOSX
    #define USE_WXMACFONTDLG 1
#else
    #define USE_WXMACFONTDLG 0
#endif

#ifdef __WXGTK__
    #define USE_WXGTK 1
#else
    #define USE_WXGTK 0
#endif

#define USE_GENERIC_DIALOGS (!USE_WXUNIVERSAL && !USE_DLL)

#define USE_COLOURDLG_GENERIC \
    ((USE_WXMSW || USE_WXMAC) && USE_GENERIC_DIALOGS && wxUSE_COLOURDLG)
#define USE_DIRDLG_GENERIC \
    ((USE_WXMSW || USE_WXMAC) && USE_GENERIC_DIALOGS && wxUSE_DIRDLG)
#define USE_FONTDLG_GENERIC \
    ((USE_WXMSW || USE_WXMACFONTDLG) && USE_GENERIC_DIALOGS && wxUSE_FONTDLG)

// Turn USE_MODAL_PRESENTATION to 0 if there is any reason for not presenting difference
// between modal and modeless dialogs (ie. not implemented it in your port yet)
#if !wxUSE_BOOKCTRL
    #define USE_MODAL_PRESENTATION 0
#else
    #define USE_MODAL_PRESENTATION 1
#endif


// Turn USE_SETTINGS_DIALOG to 0 if supported
#if wxUSE_BOOKCTRL
    #define USE_SETTINGS_DIALOG 1
#else
    #define USE_SETTINGS_DIALOG 0
#endif


// ----------------------------------------------------------------------------
// resources
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// private classes
// ----------------------------------------------------------------------------

// Define a new application type, each program should derive a class from wxApp
class PaHmmApp : public wxApp
{
public:
    // override base class virtuals
    // ----------------------------

    // this one is called on application startup and is a good place for the app
    // initialization (doing it here and not in the ctor allows to have an error
    // return: if OnInit() returns false, the application terminates)
    virtual bool OnInit() wxOVERRIDE;


};

// Define a new frame type: this is going to be our main frame
class PaHmmFrame : public wxFrame
{
public:
    // ctor(s)
    PaHmmFrame(const wxString& title);

    // event handlers (these functions should _not_ be virtual)
    void OnQuit(wxCommandEvent& event);
    void OnAbout(wxCommandEvent& event);

    void OnFileOpenSingle(wxCommandEvent& event);
    void OnFileOpenMultiple(wxCommandEvent& event);

    void OnGTRModelSelect(wxCommandEvent& event);
    void OnHKYModelSelect(wxCommandEvent& event);
    void OnLGModelSelect(wxCommandEvent& event);
    void OnGammaSelect(wxCommandEvent& event);
    void OnModelParams(wxCommandEvent& event);
    void OnIndelParams(wxCommandEvent& event);

    void FileOpen(wxCommandEvent& event);
    void FilesOpen(wxCommandEvent& event);



private:
    // any class wishing to process wxWidgets events must use this macro
    wxDECLARE_EVENT_TABLE();
};

// ----------------------------------------------------------------------------
// constants
// ----------------------------------------------------------------------------

// IDs for the controls and the menu commands
enum
{
    // menu items
    PaHmm_Quit = wxID_EXIT,

    // it is important for the id corresponding to the "About" command to have
    // this standard value as otherwise it won't be handled properly under Mac
    // (where it is special and put into the "Apple" menu)
    PaHmm_About = wxID_ABOUT,
    PaHmm_FileOpenSingle,
    PaHmm_FileOpenMultiple,
    PaHmm_GTRModelSelect,
    PaHmm_HKYModelSelect,
    PaHmm_LGModelSelect,
    PaHmm_GammaSelect,
    PaHmm_IndelParamsSelect,
    PaHmm_ModelParamsSelect
};

class MyListStoreDerivedModel : public wxDataViewListStore
{
public:
    virtual bool IsEnabledByRow(unsigned int row, unsigned int col) const wxOVERRIDE;
};
#endif
