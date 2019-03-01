"row number
set number
"tags
set tags=/myrepository/tags
"hightlight
syntax enable
syntax on
colorscheme desert
"support of 
filetype plugin on
filetype plugin indent on

au GUIEnter * simalt ~x

set ambiwidth=double
set fileencoding=gb18030
set fileencodings=utf-8,gb18030,utf-16,big5

set writebackup
set nobackup

"set list
"set listchars=tab:\|\

let g:indentLine_char='|'
let g:indentLine_enabled = 1

set guifont=Consolas:h10:cANSI
set encoding=utf8

set autochdir
set nowrap

set guioptions+=b


let fortran_have_tabs=1
set tabstop=6
set cindent shiftwidth=6
set autoindent shiftwidth=6

let fortran_fold=1
set foldmethod=syntax
set foldlevelstart=99
set foldcolumn=4
