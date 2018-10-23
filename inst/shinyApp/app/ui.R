library(shiny)
library(shinydashboard)
library(shinyjs)

ui <- dashboardPage(
    
    dashboardHeader
    (
        title="Clustering Tool"
    ),
    
    dashboardSidebar
    (
        sidebarMenu
        (
            id="tabs",
            menuItem("Files Selection", tabName="t_1"),
            menuItem("Clustering", tabName="t_2"),
            menuItem("Download", tabName = "t_3")
        )
    ),
    
    dashboardBody
    (
        useShinyjs(),
        tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")),
        tabItems
        (
            tabItem
            (
                tabName="t_1",
                h2("Files Selection"),
                fluidRow
                (
                    id="t_1_fr",
                    box
					(
						id="t_1_1",width=2,height = "27vh",
						actionButton("t_1_1_add", "Add files",style="width:90%;margin-left:4.8%"),
					    actionButton("t_1_1_rm", "Remove Selection",style="width:90%;margin-left:4.8%"),
					    actionButton("t_1_1_compensate", "Compensate Selection",style="width:90%;margin-left:4.8%"),
					    actionButton("t_1_1_select_all", "Select All",style="width:90%;margin-left:4.8%;margin-top:5.8vh"),
					    actionButton("t_1_1_deselect_all", "Deselect All",style="width:90%;margin-left:4.8%")
					),
                    box
                    (
                        id="t_1_2",width=3,height = "27vh",
                        selectInput("t_1_2_select_transform", "Select Transformation", choices=list("logicle"=1, "arcsinh"=2), selected = 1),
                        box(id="t_1_2_sel_box",width=12,style="overflow:auto;height:9vh"),
                        actionButton("t_1_2_transform", "Transform Selection",style="width:90%;margin-left:4.8%")
                    ),
                    hidden(fluidRow
                    (
                        id="t_1_4", width=6,
                        box
                        (
                            height = "25vh",style="overflow:auto",
                            tableOutput("t_1_4_fileInfo")
                        )
                    )),
                    fluidRow
                    (
                        id="t_1_3"
                    )
                )
            ),
			
            tabItem
            (
                tabName="t_2",
                h2("Clustering"),
                shinyjs::hidden(fluidRow
                (
                    id="t_2_fr",style="padding:2%;width:90%;border:solid gray;border-width:2px 0px 0px 1px",
                    box
					(
						height = "10vh", width=8, id="t_2_1", style="overflow:auto;max-height:20vh",
						selectInput("t_2_1_sel", "Select Algortihms", choices = NULL, multiple = T)
					),
                    box
                    (
                        height = "10vh", width=2, id="t_2_1", style="overflow:auto;max-height:20vh",
                        actionButton("t_2_1_run", "Run Algorithms", width="90%", style="margin-top:3%;")
                    ),
                    ##NO USE FOR NOW=============================================================================================================
					# box
					# (
					# 	height = "20vh", width=4, id="t_2_2",
					# 	h4("Select Analysis Order")
					# ),
                    ##===========================================================================================================================
                    fluidRow
                    (
                        id="t_2_3"
                    )
                ))
            ),
            
            tabItem
            (
                tabName="t_3",
                h2("Download Files"),
                fluidRow
                (
                    id="t_3_fr",
                    box
                    (
                        id="t_3_1",width=2,height = "27vh",
                        selectInput("t_3_1_select_transform", "Select Transformation", choices=list("logicle"=1, "arcsinh"=2), selected = 1),
                        box(id="t_3_1_sel_box",width=12,style="overflow:auto;height:9vh"),
                        actionButton("t_3_1_itransform", "Detransformation Selection",style="width:90%;margin-left:4.8%")
                    ),
                    box
                    (
                        id="t_3_2",width=3,height = "27vh",
                        actionButton("t_3_2_icompensate", "Decompensate Selection",style="width:80%;margin-left:9.8%;margin-top:1vh"),
                        actionButton("t_3_2_select_all", "Select All", style="width:80%;margin-left:9.8%;margin-top:5vh"),
                        actionButton("t_3_2_deselect_all", "Deselect All", style="width:80%;margin-left:9.8%"),
                        downloadButton("t_3_3_dl", "Download Enriched Files", style="width:80%;margin-left:9.8%")
                    ),
                    hidden(fluidRow
                    (
                       id="t_3_4", width=6,
                       box
                       (
                           height = "25vh",style="overflow:auto",
                           tableOutput("t_3_4_fileInfo")
                       )
                    )),
                    fluidRow
                    (
                        id="t_3_4"
                    )
                )
            )
        )
        
    )
    
)