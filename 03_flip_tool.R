rm(list = ls())

library(ape)
library(phangorn)
library(shiny)

# 设置进化树路径
tree1_path <- "promoter_seq.nwk" #cds.nwk,promoter_seq.nwk,promoter_elements_tree40.nwk
tree2_path <- "promoter_elements_tree40.nwk"

tree1 <- read.tree(tree1_path)
tree2 <- read.tree(tree2_path)

file1_name <- basename(tree1_path)
file2_name <- basename(tree2_path)

# 提取共同叶节点并把非共同节点修剪掉
common_tips <- intersect(tree1$tip.label, tree2$tip.label)
if(length(common_tips) == 0) stop("错误：两个树没有共同的叶节点！")

tree1_common <- keep.tip(tree1, common_tips)
tree2_common <- keep.tip(tree2, common_tips)

assoc_matrix <- cbind(
  tree1_common$tip.label[match(common_tips, tree1_common$tip.label)],
  tree2_common$tip.label[match(common_tips, tree2_common$tip.label)]
)

# UI界面
ui <- fluidPage(
  tags$head(tags$style(HTML("
    .node-info { background-color: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 15px; }
    .plot-container { background-color: white; border-radius: 10px; padding: 15px; box-shadow: 0 4px 8px rgba(0,0,0,0.1); margin-bottom: 20px; }
    .stat-box { background-color: #e9ecef; padding: 12px; border-radius: 6px; margin-top: 10px; }
  "))),
  titlePanel(h1("进化树交互式对齐工具", style = "font-size: 28px;")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      div(class = "node-info",
          h3("操作说明"),
          p("● 从下方下拉菜单中选择要翻转的节点编号"),
          p("● 节点编号显示在下方树中蓝色圆圈内"),
          p("● 被选中的节点会变为黄色圆圈"),
          p("● 点击【翻转选中节点】按钮进行翻转")
      ),
      hr(),
      selectInput("node1", "树1节点编号", choices = NULL, width = "100%"),
      selectInput("node2", "树2节点编号", choices = NULL, width = "100%"),
      br(),
      actionButton("flip_selected1", "🔄 翻转树1选中节点", width = "100%",
                   style = "margin-bottom: 10px; background-color: #cce5ff; font-size: 18px; padding: 12px;"),
      actionButton("flip_selected2", "🔄 翻转树2选中节点", width = "100%",
                   style = "margin-bottom: 10px; background-color: #cce5ff; font-size: 18px; padding: 12px;"),
      hr(),
      actionButton("reset", "⟲ 重置所有", width = "100%",
                   style = "margin-bottom: 10px; background-color: #fff3cd; font-size: 18px; padding: 12px;"),
      hr(),
      textInput("save_filename", "保存文件名（不含扩展名）", value = "tree_comparison"),
      actionButton("save", "💾 保存图片", width = "100%",
                   style = "margin-bottom: 5px; background-color: #d4edda; font-size: 18px; padding: 12px;"),
      hr(),
      div(class = "stat-box",
          h4("最近操作"),
          verbatimTextOutput("last_action")
      )
    ),
#三个展示界面
    mainPanel(
      width = 9,
      div(class = "plot-container",
          h3("进化树对比图"),
          plotOutput("treePlot", height = "750px", width = "100%")
      ),
      fluidRow(
        column(6, div(class = "plot-container",
                      h3("树1"),
                      plotOutput("nodePlot1", height = "600px", width = "100%"))),
        column(6, div(class = "plot-container",
                      h3("树2"),
                      plotOutput("nodePlot2", height = "600px", width = "100%")))
      )
    )
  )
)

server <- function(input, output, session) {
  
  current_tree1 <- reactiveVal(tree1_common)
  current_tree2 <- reactiveVal(tree2_common)
  selected_node1 <- reactiveVal(NULL)
  selected_node2 <- reactiveVal(NULL)
  last_action <- reactiveVal("等待操作...")
  
  # 更新下拉菜单
  observe({
    req(current_tree1())
    n <- length(current_tree1()$tip.label)
    nodes <- if(current_tree1()$Nnode > 0) (n+1):(n+current_tree1()$Nnode) else NULL
    updateSelectInput(session, "node1", choices = nodes, selected = selected_node1())
  })
  observe({
    req(current_tree2())
    n <- length(current_tree2()$tip.label)
    nodes <- if(current_tree2()$Nnode > 0) (n+1):(n+current_tree2()$Nnode) else NULL
    updateSelectInput(session, "node2", choices = nodes, selected = selected_node2())
  })
  
  observeEvent(input$node1, {
    if(!is.null(input$node1) && input$node1 != "") {
      selected_node1(as.numeric(input$node1))
      last_action(paste("✓ 选中树1节点", input$node1))
    } else {
      selected_node1(NULL)
    }
  })
  observeEvent(input$node2, {
    if(!is.null(input$node2) && input$node2 != "") {
      selected_node2(as.numeric(input$node2))
      last_action(paste("✓ 选中树2节点", input$node2))
    } else {
      selected_node2(NULL)
    }
  })
  
  observeEvent(input$flip_selected1, {
    node <- selected_node1()
    if(!is.null(node)) {
      current_tree1(rotate(current_tree1(), node))
      last_action(paste("✓ 翻转树1节点", node))
    } else {
      last_action("⚠ 请先在下拉菜单中选择树1的节点")
    }
  })
  observeEvent(input$flip_selected2, {
    node <- selected_node2()
    if(!is.null(node)) {
      current_tree2(rotate(current_tree2(), node))
      last_action(paste("✓ 翻转树2节点", node))
    } else {
      last_action("⚠ 请先在下拉菜单中选择树2的节点")
    }
  })
  
  observeEvent(input$reset, {
    current_tree1(tree1_common)
    current_tree2(tree2_common)
    selected_node1(NULL); selected_node2(NULL)
    updateSelectInput(session, "node1", selected = NULL)
    updateSelectInput(session, "node2", selected = NULL)
    last_action("✓ 已重置")
  })
  
  output$treePlot <- renderPlot({
    par(mar = c(4, 3, 5, 3))
    cophyloplot(current_tree1(), current_tree2(), assoc = assoc_matrix,
                space = 100, gap = 4, col = "darkblue", lwd = 1.8,
                use.edge.length = FALSE, show.tip.label = TRUE, font = 2, cex = 1.0)
    title(main = paste(file1_name, "vs", file2_name), line = 3, cex.main = 1.6)
  })
  
  # 树1
  output$nodePlot1 <- renderPlot({
    tree <- current_tree1()
    n_tips <- length(tree$tip.label)
    selected <- selected_node1()
    par(mar = c(4, 3, 5, 5))
    plot(tree, use.edge.length = FALSE,
         main = "树1 - 内部节点编号",
         cex = 1.0, cex.main = 1.4, font = 2, label.offset = 0.8, edge.width = 1.5)
    if(tree$Nnode > 0) {
      internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)
      normal_nodes <- setdiff(internal_nodes, selected)
      if(length(normal_nodes) > 0) {
        nodelabels(text = normal_nodes, node = normal_nodes, bg = "lightblue",
                   cex = 1.1, col = "darkred", frame = "circle", font = 2)
      }
      if(!is.null(selected) && selected %in% internal_nodes) {
        nodelabels(text = selected, node = selected, bg = "yellow",
                   cex = 1.3, col = "red", frame = "circle", font = 2)
      }
    }
    mtext("从左侧下拉菜单选择节点编号 → 点击【翻转】按钮", side = 1, line = 2, cex = 1.1, col = "blue")
    legend("bottomleft", legend = c("可选内部节点", "当前选中节点"), fill = c("lightblue", "yellow"), cex = 1.1)
  })
  
  # 树2
  output$nodePlot2 <- renderPlot({
    tree <- current_tree2()
    n_tips <- length(tree$tip.label)
    selected <- selected_node2()
    par(mar = c(4, 5, 5, 3))
    plot(tree, use.edge.length = FALSE,
         main = "树2 - 内部节点编号",
         cex = 1.0, cex.main = 1.4, font = 2, label.offset = 0.8, edge.width = 1.5)
    if(tree$Nnode > 0) {
      internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)
      normal_nodes <- setdiff(internal_nodes, selected)
      if(length(normal_nodes) > 0) {
        nodelabels(text = normal_nodes, node = normal_nodes, bg = "lightblue",
                   cex = 1.1, col = "darkred", frame = "circle", font = 2)
      }
      if(!is.null(selected) && selected %in% internal_nodes) {
        nodelabels(text = selected, node = selected, bg = "yellow",
                   cex = 1.3, col = "red", frame = "circle", font = 2)
      }
    }
    mtext("从左侧下拉菜单选择节点编号 → 点击【翻转】按钮", side = 1, line = 2, cex = 1.1, col = "blue")
    legend("bottomright", legend = c("可选内部节点", "当前选中节点"), fill = c("lightblue", "yellow"), cex = 1.1)
  })
  
  output$last_action <- renderText(last_action())
  
  # 保存图片
  observeEvent(input$save, {
    fname <- input$save_filename
    if(is.null(fname) || fname == "") fname <- "tree_comparison"
    if(!grepl("\\.pdf$", fname, ignore.case = TRUE)) fname <- paste0(fname, ".pdf")
    
    pdf(fname, width = 16, height = 9)
    par(mar = c(4, 3, 5, 3))
    cophyloplot(current_tree1(), current_tree2(), assoc = assoc_matrix,
                space = 100, gap = 4, col = "darkblue", lwd = 1.8,
                use.edge.length = FALSE, show.tip.label = TRUE, font = 2, cex = 1.0)
    title(main = paste(file1_name, "vs", file2_name), line = 3, cex.main = 1.6)
    dev.off()
    last_action(paste("✓ 图片已保存为", fname))
  })
}

runApp(list(ui = ui, server = server))