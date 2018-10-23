library(shiny)
library(shinydashboard)
library(shinyjs)
library(flowCore)
library(doSNOW)
library(parallel)


server <- function(input, output, session)
{
    useShinyjs()
    #======================================================================================================================
    #======================REACTIVE VALUES=================================================================================
    #======================================================================================================================
    
    current.project <- reactiveValues(
        fcs.files = NULL,
        fcs.files.ui.colnames = NULL,
        modified.fcs.files = NULL,
        nmb.cores = detectCores(),
        file.info.table = NULL,
        file.info.table.visible.rows = NULL
    )
    
    clustering.algorithms <- reactiveValues(
        algorithms = NULL,
        parameters = NULL
    )
    
    env.var <- reactiveValues(
        tool.wd = getwd(),
        activate.analysis = F,
        clustering.done = F
    )
    
    write.enriched.FCS <- function(fcs, fcs.path)
    {
        keywords.to.save <- names(get.keywords.with.keypart.FCS(fcs, "MAPOP_pop_label"))
        keywords.to.save <- c(unlist(keywords.to.save), names(get.keywords.with.keypart.FCS(fcs, "EXPPUR__")))
        keywords.to.save <- c(unlist(keywords.to.save), names(get.keywords.with.keypart.FCS(fcs, "RF_pop_label")))
        keywords.to.save <- c(unlist(keywords.to.save), names(get.keywords.with.keypart.FCS(fcs, "CLMETH__")))
        
        write.FCS.CIPHE(fcs,fcs.path, keywords.to.save = keywords.to.save)
    }
    
    
    
    
    
    
    
    
    
    #======================================================================================================================
    #======================================================================================================================
    #==========================================LOAD FILES==================================================================
    #======================================================================================================================
    
    update.markers.list <- function(current.section, file.id)
    {
        markers.sel <- list()
        selected.markers <- isolate(input[[paste0("t_1_3_",file.id,"_mark_sel")]])
        
        lapply(1:ncol(current.project$fcs.files[[file.id]]@exprs), function(j)
        {
            markers.sel[[current.project$fcs.files.ui.colnames[[file.id]][[j]]]] <<- j
        })
        updateSelectInput(session, paste0(current.section,"_",file.id,"_mark_sel"),label = "Markers to use",
                    choices=markers.sel,
                    selected = selected.markers)
    }
    
    observe(#LOAD FILES INFORMATION
    {
        if(length(current.project$fcs.files)>0)
        {
            lapply(1:length(current.project$fcs.files), function(f)
            {
                if(is.defined(current.project$fcs.files[[f]]) && current.project$modified.fcs.files[[f]])
                {
                    idf <- names(current.project$fcs.files)[f]
                    fcs <- current.project$fcs.files[[f]]
                    
                    #POP COL LOADING---------------------
                    pop.col.sel <- 1:ncol(fcs@exprs)
                    names(pop.col.sel) <- lapply(1:ncol(fcs@exprs), function(j)
                    {
                        d <- fcs@description[[paste0("$P",j,"S")]]
                        if(is.null(d) || !is.na(d) || d != "" || d != " ")
                        {
                            d <- current.project$fcs.files.ui.colnames[[f]][[j]]
                        }
                        names(d) <- NULL
                        
                        return(unlist(d))
                    })
                    curr.file.label <- NULL
                    if(keyword.exists.FCS(fcs,"RF_pop_label"))
                    {
                        curr.file.label <- as.numeric(get.keywords.with.keypart.FCS(fcs,"RF_pop_label")[[1]][[1]])
                    }
                    
                    #UI CREATION------------------------
                    markers.sel <- list()
                    
                    lapply(1:ncol(fcs@exprs), function(j)
                    {
                        markers.sel[[current.project$fcs.files.ui.colnames[[f]][[j]]]] <<- j
                    })
                    insertUI("#t_1_3",
                             "beforeEnd",
                             fluidRow
                             (
                                 style="padding-left:1.7vw;padding-right:0.2vw",id=paste0("t_1_3_",f,"_fr"),
                                 box
                                 (
                                     title=names(current.project$fcs.files)[f],collapsible=TRUE,width=10,collapsed=F,
                                     id=paste0("t_1_3_",f), 
                                     selectInput(paste0("t_1_3_",f,"_mark_sel"),label = "Markers to use",
                                                 choices=markers.sel,
                                                 selected=markers.sel,
                                                 multiple = T),
                                     actionButton(paste0("t_1_3_",f,"_mark_all"), "Select all"),
                                     actionButton(paste0("t_1_3_",f,"_unmark_all"), "Deselect all")
                                 ),
                                 box
                                 (
                                     width=1,height="12vh", style="padding-top:2vh",
                                     checkboxInput(paste0("t_1_3_",f,"_cbox"), "Select", value = F)
                                 )
                             )
                    )
                    insertUI("#t_3_4",
                             "beforeEnd",
                             fluidRow
                             (
                                 style="margin-left:1.7vw;padding-right:0.2vw",id=paste0("t_3_4_",f,"_fr"),
                                 box
                                 (
                                     title=names(current.project$fcs.files)[f],collapsible=TRUE,width=10,collapsed=F,
                                     id=paste0("t_3_4_",f), 
                                     selectInput(paste0("t_3_4_",f,"_mark_sel"),label = "Markers to use",
                                                         choices=markers.sel,
                                                         selected=markers.sel,
                                                         multiple = T)
                                 ),
                                 box
                                 (
                                     width=1,height="12vh", style="padding-top:2vh",
                                     checkboxInput(paste0("t_3_4_",f,"_cbox"), "Select", value = F),
                                     actionButton(paste0("t_3_4_",f,"_mark_all"), "Select all"),
                                     actionButton(paste0("t_3_4_",f,"_unmark_all"), "Deselect all")
                                 )
                             )
                    )
                    current.project$modified.fcs.files[[f]] <<- F
                    
                    
                    #SELECT ALL BUTTON---------------------------------------------
                    
                    observeEvent(input[[paste0("t_1_3_",f,"_mark_all")]],
                    {
                        markers.sel <- list()
                        lapply(1:ncol(current.project$fcs.files[[f]]@exprs), function(j)
                        {
                            markers.sel[[current.project$fcs.files.ui.colnames[[f]][[j]]]] <<- j
                        })
                        updateSelectInput(session, paste0("t_1_3_",f,"_mark_sel"),label = "Markers to use",
                                          choices = markers.sel,
                                          selected = markers.sel)
                    })
                    
                    observeEvent(input[[paste0("t_1_3_",f,"_unmark_all")]],
                    {
                        markers.sel <- list()
                        lapply(1:ncol(current.project$fcs.files[[f]]@exprs), function(j)
                        {
                            markers.sel[[current.project$fcs.files.ui.colnames[[f]][[j]]]] <<- j
                        })
                        updateSelectInput(session, paste0("t_1_3_",f,"_mark_sel"),label = "Markers to use",
                                          choices = markers.sel,
                                          selected = NULL)
                    })
                    
                    observeEvent(input[[paste0("t_3_4_",f,"_mark_all")]],
                    {
                         markers.sel <- list()
                         lapply(1:ncol(current.project$fcs.files[[f]]@exprs), function(j)
                         {
                             markers.sel[[current.project$fcs.files.ui.colnames[[f]][[j]]]] <<- j
                         })
                         updateSelectInput(session, paste0("t_3_4_",f,"_mark_sel"),label = "Markers to use",
                                           choices = markers.sel,
                                           selected = markers.sel)
                    })
                    
                    observeEvent(input[[paste0("t_3_4_",f,"_unmark_all")]],
                    {
                         markers.sel <- list()
                         lapply(1:ncol(current.project$fcs.files[[f]]@exprs), function(j)
                         {
                             markers.sel[[current.project$fcs.files.ui.colnames[[f]][[j]]]] <<- j
                         })
                         updateSelectInput(session, paste0("t_3_4_",f,"_mark_sel"),label = "Markers to use",
                                           choices = markers.sel,
                                           selected = NULL)
                    })
                }
            })
        }
    })
    
    observeEvent(input$t_1_1_add,#ADD FILES TO SESSION
    {
        shinyjs::disable("t_1_3")
        
        m <- matrix(nrow=1,ncol=2)
        m[1,1] = "FlowFrames"
        m[1,2] = "*.csv;*.fcs"
        temp.files <- choose.files(filters = m,multi = T)
        
        if(length(temp.files) > 0)
        {
            progress.bar <- Progress$new()
            progress.bar$set("LOADING FILES", value=0)
            on.exit(progress.bar$close())
            lapply(temp.files, function(f)
            {
                l <- length(f)
                x <- NULL
                if(grepl("csv",f))
                {
                    x <- as.matrix(read.csv(f))
                    x <- flowFrame(x)
                    lapply(1:ncol(x@exprs), function(i)
                    {
                        nx <- x@description[[paste0("$P",i,"S")]]
                        if(!is.null(nx) && !is.na(nx) && nx != "" && nx != " ")
                        {
                            if( is.null(current.project$fcs.files.ui.colnames) )
                            {
                                current.project$fcs.files.ui.colnames <<- list()
                            }
                            if( is.null(current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files.ui.colnames))]] ))
                            {
                                current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files.ui.colnames))]] <<- list()
                            }
                            current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files.ui.colnames))]] <<- 
                                c(current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files.ui.colnames))]],
                                  nx)
                        }
                        else
                        {
                            if( is.null(current.project$fcs.files.ui.colnames) )
                            {
                                current.project$fcs.files.ui.colnames <<- list()
                            }
                            if( is.null(current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files.ui.colnames))]] ))
                            {
                                current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files.ui.colnames))]] <<- list()
                            }
                            current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files.ui.colnames))]] <<- 
                                c(current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files.ui.colnames))]],
                                  colnames(x)[i])
                        }
                    })
                    
                }
                else
                {
                    x <- read.FCS(f,emptyValue = FALSE)
                    lapply(1:ncol(x@exprs), function(i)
                    {
                        nx <- x@description[[paste0("$P",i,"S")]]
                        if(!is.null(nx) && !is.na(nx) && nx != "" && nx != " ")
                        {
                            if( is.null(current.project$fcs.files.ui.colnames) )
                            {
                                current.project$fcs.files.ui.colnames <<- list()
                            }
                            if( is.null(current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files))]] ))
                            {
                                current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files))]] <<- list()
                            }
                            current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files))]] <<- 
                                c(current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files))]],
                                  nx)
                        }
                        else
                        {
                            if(!is.null(nx) && !is.na(nx) && nx != "" && nx != " ")
                            {
                                if( is.null(current.project$fcs.files.ui.colnames) )
                                {
                                    current.project$fcs.files.ui.colnames <<- list()
                                }
                                if( is.null(current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files))]] ))
                                {
                                    current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files))]] <<- list()
                                }
                                current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files))]] <<- 
                                    c(current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files))]],
                                      nx)
                            }
                            else
                            {
                                if( is.null(current.project$fcs.files.ui.colnames) )
                                {
                                    current.project$fcs.files.ui.colnames <<- list()
                                }
                                if( is.null(current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files))]] ))
                                {
                                    current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files))]] <<- list()
                                }
                                current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files))]] <<- 
                                    c(current.project$fcs.files.ui.colnames[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files))]],
                                      colnames(x)[i])
                            }
                        }
                    })
                }
                
                if( is.null(current.project$fcs.files) )
                {
                    current.project$fcs.files <<- list()
                    current.project$modified.fcs.files <<- list()
                }
                current.project$fcs.files[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files))]] <<- x
                current.project$modified.fcs.files[[paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files))]] <<- T
                progress.bar$inc(1/length(temp.files),detail=paste0("adding ",f))
                
                file.vec <- matrix(ncol=3,nrow=1)
                file.vec[1,1] <- paste0(basename(substr(f,1,nchar(f)-4)),"_",length(current.project$fcs.files)-1)
                file.vec[1,2] <- trunc(file.size(f)/1024/1024*1000)/1000
                file.vec[1,3] <- ncol(x)
                
                current.project$file.info.table <<- rbind(current.project$file.info.table,
                                                          file.vec)
                current.project$file.info.table.visible.rows <<- c(current.project$file.info.table.visible.rows,
                                                                   T)
                colnames(current.project$file.info.table) <<- c("Filename", "Size (Mo)", "Number of markers")
            })
            progress.bar$close()
        }
        
        shinyjs::enable("t_1_3")
    })
    
    observeEvent(input$t_1_1_rm,#REMOVE SELECTED FILES
    {
        progress.bar <- Progress$new()
        progress.bar$set("REMOVING FILES", value=0)
        on.exit(progress.bar$close())
        if( length(current.project$fcs.files) >0 )
        {
            nmb.files.to.remove <- 0
            lapply(1:length(current.project$fcs.files), function(i)
            {
                if(input[[paste0("t_1_3_",i,"_cbox")]])
                {
                    nmb.files.to.remove <<- nmb.files.to.remove+1
                    
                }
            })
            lapply(1:length(current.project$fcs.files), function(i)
            {
                if(input[[paste0("t_1_3_",i,"_cbox")]])
                {
                    current.project$fcs.files[[i]] <<- NA
                    current.project$modified.fcs.files[[i]] <<- F
                    current.project$fcs.files.ui.colnames[[i]] <<- NA
                    progress.bar$inc(1/nmb.files.to.remove, detail=paste0("File ", i, " removed"))
                    removeUI(paste0("#t_3_4_",i,"_fr"))
                    removeUI(paste0("#t_1_3_",i,"_fr"))
                    current.project$file.info.table.visible.rows[[i]] <<-  F
                }
            })
        }
    })
    
    observeEvent(input$t_1_1_compensate,#COMPENSATE SELECTED FILES
    {
         selected.files <- 0
         if( length(current.project$fcs.files) >0 )
         {
            lapply(1:length(current.project$fcs.files), function(i)
            {
                if(input[[paste0("t_1_3_",i,"_cbox")]] && is.defined(current.project$fcs.files[[i]]))
                {
                    selected.files <<- selected.files + 1
                }
            })
         }
         
         if(selected.files>0)
         {
             progress.bar <- Progress$new()
             progress.bar$set("COMPENSATION", value=0)
             progress.bar$inc(0,detail="please wait")
             if( length(current.project$fcs.files) >0 )
             {
                 tmp.fcs.files <- isolate(reactiveValuesToList(current.project))$fcs.files
                 tmp.fcs.files.names <- names(isolate(reactiveValuesToList(current.project))$fcs.files)
                 tmp.input <- isolate(reactiveValuesToList(input))
                 
                 files.sizes <- unlist(sapply(tmp.fcs.files, function(curr.f){return(object.size(curr.f))}))
                 nmb.cl <- get.nmb.cores.max(files.sizes, available.cores = current.project$nmb.cores, x.cores = 0.1,
                                             x.ram = 0.3, correction.coef = 1.05, separate.by.files = T)
                 cl <- makeCluster(nmb.cl)
                 registerDoSNOW(cl)
                 
                 progress.fct <- function(i)
                 {
                     par.name <- names(tmp.fcs.files)[i]
                     progress.bar$inc(1/selected.files,
                                      detail=paste0(par.name))
                 }
                 
                 in.time <- Sys.time()
                 tmp.fcs.files <- foreach(f.id=1:length(tmp.fcs.files), 
                                          .options.snow = list(progress=progress.fct), 
                                          .packages = c("flowCore"),
                                          .export = c("m.compensate","is.defined")) %dopar%
                 {
                     fcs <- tmp.fcs.files[[f.id]]
                     if(is.defined(fcs))
                     {
                         if(tmp.input[[paste0("t_1_3_",f.id,"_cbox")]])
                         {
                             fcs <- m.compensate(fcs)
                         }
                     }
                     return(fcs)
                 }
                 print("EXEC TIME: ")
                 print(Sys.time()-in.time)
                 
                 stopCluster(cl)
                 current.project$fcs.files <<- tmp.fcs.files
                 names(current.project$fcs.files) <<- tmp.fcs.files.names
                 
             }
             progress.bar$close()
         }
         else
         {
             progress.bar <- Progress$new()
             progress.bar$set("NOTHING TO BE DONE", value=1)
             delay(1500, progress.bar$close())
         }
         
    })
    
    observeEvent(input$t_1_2_transform,#TRANSFORM SELECTED FILES
    {
         selected.files <- 0
         if( length(current.project$fcs.files) >0 )
         {
             lapply(1:length(current.project$fcs.files), function(i)
             {
                 if(input[[paste0("t_1_3_",i,"_cbox")]] && is.defined(current.project$fcs.files[[i]]))
                 {
                     selected.files <<- selected.files+1
                 }
             })
         }
         
         if(selected.files>0)
         {
             progress.bar <- Progress$new()
             progress.bar$set("TRANSFORMATION: ", value=0)
             progress.bar$inc(0,detail="please wait")
             if( length(current.project$fcs.files) >0 )
             {
                 tmp.fcs.files <- isolate(reactiveValuesToList(current.project))$fcs.files
                 tmp.fcs.files.names <- names(isolate(reactiveValuesToList(current.project))$fcs.files)
                 tmp.input <- isolate(reactiveValuesToList(input))
                 
                 selected.algo <- m.transform.logicle
                 selected.algo.params <- NULL
                 if(is.defined(input$t_1_2_select_transform) && input$t_1_2_select_transform != "" && 
                    input$t_1_2_select_transform != " ")
                 {
                     selected.transform <- as.numeric(input$t_1_2_select_transform)
                     if(selected.transform == 1)
                     {
                         selected.algo <- m.transform.logicle
                         selected.algo.params <- NULL
                     }
                     else
                     {
                         selected.algo <- m.transform.asinh
                         selected.algo.params <- as.numeric(tmp.input$t_1_2_sel_arcsinh)
                     }
                 }
                 
                 files.sizes <- unlist(sapply(tmp.fcs.files, function(curr.f){return(object.size(curr.f))}))
                 nmb.cl <- get.nmb.cores.max(files.sizes, available.cores = current.project$nmb.cores, x.cores = 0.1,
                                             x.ram = 0.3, correction.coef = 1.05, separate.by.files = T)
                 cl <- makeCluster(nmb.cl)
                 registerDoSNOW(cl)
                 
                 progress.fct <- function(i)
                 {
                     par.name <- names(tmp.fcs.files)[i]
                     progress.bar$inc(1/selected.files,
                                      detail=paste0(par.name))
                 }
                 
                 in.time <- Sys.time()
                 tmp.fcs.files <- foreach(f.id=1:length(tmp.fcs.files), 
                                          .options.snow = list(progress=progress.fct), 
                                          .packages = c("flowCore"),
                                          .export = c("selected.algo", "selected.algo.params", "is.defined")) %dopar%
                 {
                      fcs <- tmp.fcs.files[[f.id]]
                      if(is.defined(fcs))
                      {
                          if(tmp.input[[paste0("t_1_3_",f.id,"_cbox")]])
                          {
                              fcs.col <- colnames(fcs@exprs)[as.numeric(tmp.input[[paste0("t_1_3_",f.id,"_mark_sel")]])]
                              fcs <- selected.algo(fcs, fcs.col, selected.algo.params)
                          }
                      }
                      return(fcs)
                 }
                 print("EXEC TIME: ")
                 print(Sys.time()-in.time)
                 
                 stopCluster(cl)
                 current.project$fcs.files <<- tmp.fcs.files
                 names(current.project$fcs.files) <<- tmp.fcs.files.names
                 
             }
             progress.bar$close()
         }
         else
         {
             progress.bar <- Progress$new()
             progress.bar$set("NOTHING TO BE DONE", value=1)
             delay(1500, progress.bar$close())
         }
         
    })
    
    observe(#CHANGE TRANSFORM UI
    {
        if(is.defined(input$t_1_2_select_transform) && input$t_1_2_select_transform != "" && 
           input$t_1_2_select_transform != " ")
        {
            selected.transform <- as.numeric(input$t_1_2_select_transform)
            if(selected.transform == 1)
            {
                removeUI("#t_1_2_sel_fr")
            }
            else
            {
                insertUI("#t_1_2_sel_box",
                         "beforeEnd",
                         fluidRow
                         (
                             id="t_1_2_sel_fr", style="width:90%;margin-left:4.8%",
                             textInput("t_1_2_sel_arcsinh", "Arcsinh cofactor", value = "5")
                         )
                )
            }
        }
    })
    
    observe(#ACTIVATE UI
    {
        if(env.var$activate.analysis)
        {
            shinyjs::enable("t_1_1_rm")
            shinyjs::enable("t_1_1_compensate")
            shinyjs::enable("t_1_1_select_all")
            shinyjs::enable("t_1_1_deselect_all")
            shinyjs::enable("t_1_2_transform")
        }
        else
        {
            shinyjs::disable("t_1_1_rm")
            shinyjs::disable("t_1_1_compensate")
            shinyjs::disable("t_1_1_select_all")
            shinyjs::disable("t_1_1_deselect_all") 
            shinyjs::disable("t_1_2_transform")
        }
    })
    
    observeEvent(input$t_1_1_select_all,#SELECT ALL FILES
    {
        if( length(current.project$fcs.files) >0 )
        {
            lapply(1:length(current.project$fcs.files), function(f)
            {
                if(is.defined(current.project$fcs.files[[f]]))
                {
                    updateCheckboxInput(session, paste0("t_1_3_",f,"_cbox"), value = T)
                }
            })
        }
    })
    
    observeEvent(input$t_1_1_deselect_all,#DESELECT ALL FILES
    {
        if( length(current.project$fcs.files) >0 )
        {
            lapply(1:length(current.project$fcs.files), function(f)
            {
                if(is.defined(current.project$fcs.files[[f]]))
                {
                    updateCheckboxInput(session, paste0("t_1_3_",f,"_cbox"), value = F)
                }
            })
        }
    })
    
    file.table.fct <- function()
    {
        l <- sum(unlist(current.project$file.info.table.visible.rows))
        tmp.mat <- current.project$file.info.table[unlist(current.project$file.info.table.visible.rows),]
        if(l==1)
        {
            tmp.mat <- t(tmp.mat)
            colnames(current.project$file.info.table) <<- c("Filename", "Size (Mo)", "Number of markers")
        }
        
        return(tmp.mat)
        
    }
    
    output$t_1_4_fileInfo <- renderTable(file.table.fct()) #FILES INFORMATION
    output$t_3_4_fileInfo <- renderTable(file.table.fct()) #FILES INFORMATION
    
    observe(#SHOW/HIDE FILES INFORMATION
    {
        if(!is.null(current.project$file.info.table))
        {
            l <- sum(unlist(current.project$file.info.table.visible.rows))
            if(l > 0)
            {
                
                shinyjs::show("t_3_4")
                shinyjs::show("t_1_4")
            }
            else
            {
                shinyjs::hide("t_3_4")
                shinyjs::hide("t_1_4")
            }
        }
        else
        {
            shinyjs::hide("t_1_4")
            shinyjs::hide("t_3_4")
        }
    })
    
    
    
    
    
    
    
    
    #======================================================================================================================
    #======================================================================================================================
    #==========================================CLUSTERING==================================================================
    #======================================================================================================================
    
    observe(#ACTIVATE BUTTON
    {
        if(length(current.project$fcs.files)>0)
        {
            if(sum(is.na(unlist(current.project$fcs.files)))<length(current.project$fcs.files))
            {
                env.var$activate.analysis = T
            }
            else
            {
                env.var$activate.analysis = F
            }
        }
        else
        {
            env.var$activate.analysis = F
        }
    })

    observe(#SHOW UI
    {
        if(env.var$activate.analysis)
        {
            shinyjs::show("t_2_fr")
        }
        else
        {
            shinyjs::hide("t_2_fr")
        }
    })

    observe(#LISTS AVAILABLE ALGORITHMS
    {
        if( length(clustering.algorithms$algorithms)==0 )
        {
            temp.dir <- paste0(env.var$tool.wd,"/MethodsFolder/")
            methods.files <- list.files(temp.dir, pattern = ".R", full.names = F)
            if( is.null(clustering.algorithms$algorithls) )
            {
                clustering.algorithms$algorithms <- list()
            }
            if( is.null(clustering.algorithms$parameters) )
            {
                clustering.algorithms$parameters <- list()
            }

            lapply(methods.files, function(f)
            {
                source(paste0(temp.dir,f))

                clustering.algorithms$algorithms[[strsplit(f,".R", fixed = T)[[1]][1]]] <<- strsplit(f,".R", fixed = T)[[1]][1]
                clustering.algorithms$parameters[[strsplit(f,".R", fixed = T)[[1]][1]]] <<- fct.parameters
            })
        }

        if(length(clustering.algorithms$algorithms)>0)
        {
            algo <- names(clustering.algorithms$algorithms)
            names(algo) <- algo
            updateSelectInput(session, "t_2_1_sel", "Select Algorithms", choices=algo, selected = algo)
        }
    })

    observe(#PLOTS UI DEPENDING ON SELECTED ALGORITHMS
    {
        removeUI("#t_2_3")
        if( length(current.project$fcs.files)>0 && is.defined(input$t_2_1_sel) )
        {
            insertUI("#t_2_fr",
                     "beforeEnd",
                     fluidRow
                     (
                         id="t_2_3", style="padding-left:1.2%"
                     )
            )

            lapply(1:length(input$t_2_1_sel), function(k)
            {
                insertUI("#t_2_3",
                         "beforeEnd",
                         box
                         (
                             title=paste0(input$t_2_1_sel[[k]], ": Parameters"), id=paste0("t_2_3_",k),style="padding:2vw",
                             width=10, collapsible=T, collapsed=F
                         )
                )

                if( !is.null(clustering.algorithms$parameters[[input$t_2_1_sel[[k]]]]) )
                {
                    lapply(1:length(clustering.algorithms$parameters[[input$t_2_1_sel[[k]]]]), function(p)
                    {
                        par <- clustering.algorithms$parameters[[input$t_2_1_sel[[k]]]][[p]]
                        par.name <- names(clustering.algorithms$parameters[[input$t_2_1_sel[[k]]]])[p]

                        insertUI(paste0("#t_2_3_",k),
                                 "beforeEnd",
                                 sliderInput(paste0("t_2_3_",k,"_",p),par.name,min = as.numeric(par[1]),max=as.numeric(par[3]),
                                             step=as.numeric(par[2]),value=c(as.numeric(par[4]),as.numeric(par[4])))
                        )
                    })
                }

            })

        }
    })

    #NO USE FOR NOW========================================================================================================

    # observe(#UPDATES ANALYSIS ORDER UI
    # {
    #     if( length(current.project$fcs.files)>0 && is.defined(input$t_2_1_sel) )
    #     {
    #         removeUI("#t_2_2_dropBox")
    #         insertUI("#t_2_2",
    #                  "beforeEnd",
    #                  dropUI("t_2_2_dropBox", style = "height:100%", col_n = 3)
    #         )
    #         lapply(1:length(input$t_2_1_sel), function(alg.ID)
    #         {
    #             insertUI("#t_2_2_dropBox",
    #                      "beforeEnd",
    #                      dragUI(paste0("t_2_2_",alg.ID,"drag"),
    #                             h4(input$t_2_1_sel[[alg.ID]]),
    #                             style="width:90%")
    #             )
    #         })
    #     }
    # })

    # observeEvent(input[["t_2_2_dropBox"]],#CHECK THE ELEMENTS IN THE ANALYSIS ORDER DROPBOX (temporary)
    # {
    #     y <- input[["t_2_2_dropBox"]]
    #     print(y)
    # })

    #=====================================================================================================================

    observeEvent(input$t_2_1_run,#RUNS ANALYSES
    {
        if(length(current.project$fcs.files)>0)
        {
            if(is.defined(input$t_2_1_sel) && length(input$t_2_1_sel)>0)
            {
                lapply(1:length(input$t_2_1_sel), function(alg.id)
                {
                    curr.algo <- input$t_2_1_sel[[alg.id]]
                    if(length(clustering.algorithms$parameters[[curr.algo]]) > 0)
                    {
                        params <- list()
                        params <- lapply(1:length(clustering.algorithms$parameters[[curr.algo]]), function(p)
                        {
                            x <- list(1,2,3)
                            values <- c()
                            if(is.defined(input[[paste0("t_2_3_",alg.id,"_",p)]]))
                            {
                                tmp <- input[[paste0("t_2_3_",alg.id,"_",p)]]
                                x[[1]] <- as.numeric(tmp)[[1]]
                                x[[2]] <- as.numeric(tmp)[[2]]
                                x[[3]] <- as.numeric(clustering.algorithms$parameters[[curr.algo]][[p]][[2]])
                            }
                            if(is.defined(x[[3]]))
                            {
                                values <- seq(as.numeric(x[[1]]),as.numeric(x[[2]]),as.numeric(x[[3]]))
                            }
                            return(values)
                        })
                        names(params) <- names(clustering.algorithms$parameters[[curr.algo]])
                    }
                    runs.params.list <- run.algo.param.combi(params)

                    #DEFINES PARALLEL WORK=========================================================================
                    nmb.runs <- length(runs.params.list)
                    runs.per.core <- nmb.runs/current.project$nmb.cores

                    #Store reactive and global values
                    tmp.curr.proj <- isolate(reactiveValuesToList(current.project))
                    tmp.input <- isolate(reactiveValuesToList(input))
                    tmp.tool.wd <- isolate(reactiveValuesToList(env.var))$tool.wd
                    tmp.algo.params <- isolate(reactiveValuesToList(clustering.algorithms))$parameters
                    
                    tmp.L1 <- lapply(1:length(runs.params.list), function(tmp.id)
                    {
                        return(list(runs.params.list[[tmp.id]], tmp.id))
                    })
                    
                    tmp.L2 <- lapply(1:length(tmp.curr.proj$fcs.files), function(tmp.id)
                    {
                        return(list(tmp.curr.proj$fcs.files[[tmp.id]], names(tmp.curr.proj$fcs.files)[tmp.id], tmp.id))
                    })
                    
                    tmp.foreach.list <- run.algo.combi(list(tmp.L1, tmp.L2))
                    L1.1 <- lapply(1:length(tmp.foreach.list), function(i)
                    {
                        return(tmp.foreach.list[[i]][[1]])
                    })
                    L1.2 <- lapply(1:length(tmp.foreach.list), function(i)
                    {
                        return(tmp.foreach.list[[i]][[2]])
                    })
                    L2.1 <- lapply(1:length(tmp.foreach.list), function(i)
                    {
                        return(tmp.foreach.list[[i]][[3]])
                    })
                    L2.2 <- lapply(1:length(tmp.foreach.list), function(i)
                    {
                        return(tmp.foreach.list[[i]][[4]])
                    })
                    L2.3 <- lapply(1:length(tmp.foreach.list), function(i)
                    {
                        return(tmp.foreach.list[[i]][[5]])
                    })
                    
                    files.sizes <- unlist(sapply(tmp.curr.proj$fcs.files, function(curr.f){return(object.size(curr.f))}))
                    nmb.cl <- get.nmb.cores.max(files.sizes, available.cores = current.project$nmb.cores, x.cores = 0.1,
                                                x.ram = 0.3, correction.coef = 1.05, separate.by.files = F)
                    cl <- makeCluster(nmb.cl)
                    registerDoSNOW(cl)

                    progress.bar <- Progress$new()
                    progress.bar$set(paste0("ALGORITHM: ", curr.algo), value=0)
                    progress.bar$inc(0,detail="fetching parameters")
                    progress.bar.fct <- function(i)
                    {
                        par.val <- unlist(L1.1[[as.integer(i)]])
                        par.name <- names(tmp.algo.params[[curr.algo]])
                        progress.bar$inc(1/length(L1.1),
                                         detail=paste0(par.name,": ",par.val))
                    }
                    in.time <- Sys.time()
                    
                    temp.out <- foreach(run.id=L1.2, run.parameters.values=L1.1, fcs=L2.1, fcs.name=L2.2, f.id=L2.3,
                                        .options.snow = list(progress=progress.bar.fct),
                                        .packages=c("flowCore","microbenchmark"),
                                        .export = c("is.defined","benchmark.method","benchmark.source.method","add.keyword.to.fcs","alg.id",
                                                    "curr.algo","enrich.FCS", "params","tmp.input","tmp.tool.wd","tmp.algo.params")) %dopar%
                    {
                        added.keyword <- NULL
                        added.keyword.name <- NULL
                        if(is.defined(fcs))
                        {
                            markers_col <- 1
                            if( is.defined(tmp.input[[paste0("t_1_3_",f.id,"_mark_sel")]]) )
                            {
                                markers_col <- tmp.input[[paste0("t_1_3_",f.id,"_mark_sel")]]
                            }
                            
                            benchmark.source.method(paste0(tmp.tool.wd,"/MethodsFolder/"),curr.algo)
                            method.output <- benchmark.method(curr.algo, fcs, run.parameters.values, markers_col)
                            
                            tmp.labels <- method.output[[1]]
                            fcs <- enrich.FCS(fcs,tmp.labels)
                            
                            added.keyword <- paste0("CLMETH__",curr.algo,"__",ncol(fcs@exprs),"__")
                            markers.txt <- "NULL"
                            if(length(markers_col) > 0)
                            {
                                markers.txt <- ""
                                lapply(1:length(markers_col), function(i)
                                {
                                    markers.txt <<- paste0(markers.txt,markers_col[[i]],".-.")
                                })
                            }
                            params.txt <- "NULL"
                            if(length(run.parameters.values) > 0)
                            {
                                params.txt <- ""
                                lapply(1:length(run.parameters.values), function(i)
                                {
                                    p.name <- names(params)[i]
                                    p.val <- run.parameters.values[[i]]
                                    params.txt <<- paste0(params.txt,p.name,"-",p.val,".-.")
                                })
                            }
                            added.keyword <- paste0(added.keyword,markers.txt,"__",params.txt)
                            added.keyword.name <- paste0("CLMETH__",curr.algo,"__",ncol(fcs@exprs))
                            
                            return(list(fcs@exprs[,ncol(fcs@exprs)], colnames(fcs@exprs)[ncol(fcs@exprs)], 
                                        added.keyword, added.keyword.name,
                                        fcs.name))
                        }
                        else
                        {
                            return(NULL)
                        }
                    }
                    
                    progress.bar$inc(0,detail="finalizing clustering")
                    stopCluster(cl)
                    print(paste("EXECUTION TIME:"))
                    print(Sys.time()-in.time)
                    
                    fcs.files <- current.project$fcs.files
                    
                    lapply(1:length(temp.out), function(file.id)
                    {
                        if(is.defined(temp.out[[file.id]]))
                        {
                            tmp.file.labels <- matrix(temp.out[[file.id]][[1]],ncol=1)
                            
                            tmp.name <- ""
                            tmp.val <- strsplit(temp.out[[file.id]][[2]], ".", T)[[1]]
                            if(length(tmp.val)>1)
                            {
                                lapply(1:(length(tmp.val) - 1), function(s.id)
                                {
                                    tmp.name <<- paste0(tmp.name, tmp.val[[s.id]], ".")
                                })
                            }
                            tmp.name <- paste0(tmp.name, ncol(fcs.files[[ temp.out[[file.id]][[5]] ]])+1)
                            colnames(tmp.file.labels) <- tmp.name
                            
                            fcs.files[[ temp.out[[file.id]][[5]] ]] <<- enrich.FCS(fcs.files[[ temp.out[[file.id]][[5]] ]], tmp.file.labels)
                            
                            new.key <- modify.keyword.value(temp.out[[file.id]][[3]], "__", 3,
                                                            ncol(fcs.files[[ temp.out[[file.id]][[5]] ]]@exprs))
                            new.key.name <- modify.keyword.value(temp.out[[file.id]][[4]], "__", 3, 
                                                                 ncol(fcs.files[[ temp.out[[file.id]][[5]] ]]@exprs))
                            
                            fcs.files[[ temp.out[[file.id]][[5]] ]] <<- add.keyword.to.fcs(fcs.files[[ temp.out[[file.id]][[5]] ]], 
                                                                                                 new.key, 
                                                                                                 new.key.name)
                            current.project$fcs.files.ui.colnames[[ temp.out[[file.id]][[5]] ]] <<- 
                                c(current.project$fcs.files.ui.colnames[[ temp.out[[file.id]][[5]] ]],
                                  tmp.name)
                        }
                    })
                    
                    current.project$fcs.files <<- fcs.files
                    for(f in 1:length(current.project$fcs.files))
                    {
                        if(is.defined(current.project$fcs.files[[f]]))
                        {
                            update.markers.list(current.section="t_1_3",f)
                            update.markers.list(current.section="t_3_4",f)
                        }
                    }
                    
                    progress.bar$close()
                })
            }
        }
    })
    
    
    
    
    
    
    #======================================================================================================================
    #======================================================================================================================
    #==========================================DOWNLOAD FILES==============================================================
    #======================================================================================================================
    
    observeEvent(input$t_3_2_icompensate,#INVERSE COMPENSATION
    {
         selected.files <- 0
         if( length(current.project$fcs.files) >0 )
         {
             lapply(1:length(current.project$fcs.files), function(i)
             {
                 if(input[[paste0("t_3_4_",i,"_cbox")]] && is.defined(current.project$fcs.files[[i]]))
                 {
                     selected.files <<- selected.files + 1
                 }
             })
         }
         
         if(selected.files>0)
         {
             progress.bar <- Progress$new()
             progress.bar$set("INVERTING COMPENSATION", value=0)
             progress.bar$inc(0,detail="please wait")
             if( length(current.project$fcs.files) >0 )
             {
                 tmp.fcs.files <- isolate(reactiveValuesToList(current.project))$fcs.files
                 tmp.fcs.files.names <- names(isolate(reactiveValuesToList(current.project))$fcs.files)
                 tmp.input <- isolate(reactiveValuesToList(input))
                 
                 files.sizes <- unlist(sapply(tmp.fcs.files, function(curr.f){return(object.size(curr.f))}))
                 nmb.cl <- get.nmb.cores.max(files.sizes, available.cores = current.project$nmb.cores, x.cores = 0.1,
                                             x.ram = 0.3, correction.coef = 1.05, separate.by.files = T)
                 cl <- makeCluster(nmb.cl)
                 registerDoSNOW(cl)
                 
                 progress.fct <- function(i)
                 {
                     par.name <- names(tmp.fcs.files)[i]
                     progress.bar$inc(1/selected.files,
                                      detail=paste0(par.name))
                 }
                 
                 in.time <- Sys.time()
                 tmp.fcs.files <- foreach(f.id=1:length(tmp.fcs.files), 
                                          .options.snow = list(progress=progress.fct), 
                                          .packages = c("flowCore"),
                                          .export = c("m.inv.compensate","is.defined")) %dopar%
                  {
                      fcs <- tmp.fcs.files[[f.id]]
                      if(is.defined(fcs))
                      {
                          if(tmp.input[[paste0("t_3_4_",f.id,"_cbox")]])
                          {
                              fcs <- m.inv.compensate(fcs)
                          }
                      }
                      return(fcs)
                  }
                 print("EXEC TIME: ")
                 print(Sys.time()-in.time)
                 
                 stopCluster(cl)
                 current.project$fcs.files <<- tmp.fcs.files
                 names(current.project$fcs.files) <<- tmp.fcs.files.names
                 
             }
             progress.bar$close()
         }
         else
         {
             progress.bar <- Progress$new()
             progress.bar$set("NOTHING TO BE DONE", value=1)
             delay(1500, progress.bar$close())
         }
         
     })

    observeEvent(input$t_3_1_itransform,#INVERSE TRANSFORMATION
    {
         selected.files <- 0
         if( length(current.project$fcs.files) >0 )
         {
             lapply(1:length(current.project$fcs.files), function(i)
             {
                 if(input[[paste0("t_3_4_",i,"_cbox")]] && is.defined(current.project$fcs.files[[i]]))
                 {
                     selected.files <<- selected.files+1
                 }
             })
         }
         
         if(selected.files>0)
         {
             progress.bar <- Progress$new()
             progress.bar$set("INVERTING TRANSFORMATION: ", value=0)
             progress.bar$inc(0,detail="please wait")
             if( length(current.project$fcs.files) >0 )
             {
                 tmp.fcs.files <- isolate(reactiveValuesToList(current.project))$fcs.files
                 tmp.fcs.files.names <- names(isolate(reactiveValuesToList(current.project))$fcs.files)
                 tmp.input <- isolate(reactiveValuesToList(input))
                 
                 selected.algo <- m.inv.transform.logicle
                 selected.algo.params <- NULL
                 if(is.defined(input$t_3_1_select_transform) && input$t_3_1_select_transform != "" && 
                    input$t_3_1_select_transform != " ")
                 {
                     selected.transform <- as.numeric(input$t_3_1_select_transform)
                     if(selected.transform == 1)
                     {
                         selected.algo <- m.inv.transform.logicle
                         selected.algo.params <- NULL
                     }
                     else
                     {
                         selected.algo <- m.inv.transform.asinh
                         selected.algo.params <- as.numeric(tmp.input$t_3_1_sel_arcsinh)
                     }
                 }
                 
                 files.sizes <- unlist(sapply(tmp.fcs.files, function(curr.f){return(object.size(curr.f))}))
                 nmb.cl <- get.nmb.cores.max(files.sizes, available.cores = current.project$nmb.cores, x.cores = 0.1,
                                             x.ram = 0.3, correction.coef = 1.05, separate.by.files = T)
                 cl <- makeCluster(nmb.cl)
                 registerDoSNOW(cl)
                 
                 progress.fct <- function(i)
                 {
                     par.name <- names(tmp.fcs.files)[i]
                     progress.bar$inc(1/selected.files,
                                      detail=paste0(par.name))
                 }
                 
                 in.time <- Sys.time()
                 tmp.fcs.files <- foreach(f.id=1:length(tmp.fcs.files), 
                                          .options.snow = list(progress=progress.fct), 
                                          .packages = c("flowCore"),
                                          .export = c("selected.algo", "selected.algo.params", "is.defined")) %dopar%
                 {
                      fcs <- tmp.fcs.files[[f.id]]
                      if(is.defined(fcs))
                      {
                          if(tmp.input[[paste0("t_3_4_",f.id,"_cbox")]])
                          {
                              fcs.col <- colnames(fcs@exprs)[as.numeric(tmp.input[[paste0("t_3_4_",f.id,"_mark_sel")]])]
                              fcs <- selected.algo(fcs, fcs.col, selected.algo.params)
                          }
                      }
                      return(fcs)
                 }
                 print("EXEC TIME: ")
                 print(Sys.time()-in.time)
                 
                 stopCluster(cl)
                 current.project$fcs.files <<- tmp.fcs.files
                 names(current.project$fcs.files) <<- tmp.fcs.files.names
                 
             }
             progress.bar$close()
         }
         else
         {
             progress.bar <- Progress$new()
             progress.bar$set("NOTHING TO BE DONE", value=1)
             delay(1500, progress.bar$close())
         }
         
    })
    
    observeEvent(input$t_3_2_select_all,#SELECT ALL FILES
    {
         if( length(current.project$fcs.files) >0 )
         {
             lapply(1:length(current.project$fcs.files), function(f)
             {
                 if(is.defined(current.project$fcs.files[[f]]))
                 {
                     updateCheckboxInput(session, paste0("t_3_4_",f,"_cbox"), value = T)
                 }
             })
         }
    })
    
    observeEvent(input$t_3_2_deselect_all,#DESELECT ALL FILES
    {
         if( length(current.project$fcs.files) >0 )
         {
             lapply(1:length(current.project$fcs.files), function(f)
             {
                 if(is.defined(current.project$fcs.files[[f]]))
                 {
                     updateCheckboxInput(session, paste0("t_3_4_",f,"_cbox"), value = F)
                 }
             })
         }
    })
    
    observe(#CHANGE INV TRANSFORM UI
    {
        if(is.defined(input$t_3_1_select_transform) && input$t_3_1_select_transform != "" && 
           input$t_3_1_select_transform != " ")
        {
            selected.transform <- as.numeric(input$t_3_1_select_transform)
            if(selected.transform == 1)
            {
                removeUI("#t_3_1_sel_fr")
            }
            else
            {
                insertUI("#t_3_1_sel_box",
                         "beforeEnd",
                         fluidRow
                         (
                             id="t_3_1_sel_fr", style="width:90%;margin-left:4.8%",
                             textInput("t_3_1_sel_arcsinh", "Arcsinh cofactor", value = "5")
                         )
                )
            }
        }
    })
    
    observe(#ACTIVATE UI
    {
        if(env.var$activate.analysis)
        {
            shinyjs::enable("t_3_2_select_all")
            shinyjs::enable("t_3_2_deselect_all")
            shinyjs::enable("t_3_2_icompensate")
            shinyjs::enable("t_3_1_itransform")
            shinyjs::enable("t_3_3_dl")
        }
        else
        {
            shinyjs::disable("t_3_1_itransform")
            shinyjs::disable("t_3_2_select_all")
            shinyjs::disable("t_3_2_deselect_all")
            shinyjs::disable("t_3_2_icompensate")
            shinyjs::disable("t_3_3_dl")
        }
    })

    output$t_3_3_dl <- downloadHandler(
        filename = function()
        {
            paste0("output.zip")
        },
        content = function(file)
        {
            f.names <- c()
            if(length(current.project$fcs.files)>0)
            {
                lapply(1:length(current.project$fcs.files), function(f)
                {
                    if(is.defined(current.project$fcs.files[[f]]))
                    {
                        idf <- names(current.project$fcs.files)[f]
                        fcs <- current.project$fcs.files[[f]]

                        save.name <- paste0(idf,".fcs")
                        write.enriched.FCS(fcs,save.name)
                        f.names <<- c(f.names, save.name)
                    }
                })
            }
            zip(file,f.names)
            file.remove(f.names)
        }
    )
    
}

