library(tidyverse)
library(lubridate)
library(vroom)
library(rlist) # for list to tibble combined by rows
library(deSolve)
library(Cairo)
library(scales) # for the comma formatting
ggplot2::theme_set(theme_classic())

# 1. making and wrangling data --------------------------------------------
jeju <- tibble(
    date = seq(ymd("2021-05-31"), ymd("2021-06-20"), by = "days")
    # case = c(?)
)
jeju2 <- jeju %>% 
    mutate(case_MA = slider::slide_dbl(case, mean, .before = 3, .after = 3, .complete = TRUE))
init_i <- jeju2 %>%
    filter(!is.na(case_MA)) %>% 
    tail(1) %>% 
    pull(case_MA) # 7월 이동평균의 마지막값을 기준으로 infected 수를 잡고 적합 및 plotting 예정

# 2. modeling with SEIR models --------------------------------------------
## (1)  making a SEIR model function 
jeju_seir <- function(time, state, parameters){
    with(as.list(c(state, parameters)),{
        N = S1+S2+E+I1+I2+R
        lambda = beta*(I/N)
        dS1 = -lambda*S1
        dS2 = -lambda*S2
        dE = lambda*S1 + lambda*S2 - sigma*E
        dI = sigma*E - gamma*I
        dR = gamma*I
        
        return(list(c(dS1, dS2, dE, dI, dR)))
    }
    )
}

## (2) Model inputs
param_1.3 <- c(gamma = 1/30, beta = 1.3/30, sigma = 1/14) # R = 1.3
param_1.2 <- c(gamma = 1/30, beta = 1.2/30, sigma = 1/14) # R = 1.2
param_1.1 <- c(gamma = 1/30, beta = 1.1/30, sigma = 1/14) # R = 1.1
param_1.0 <- c(gamma = 1/30, beta = 1.0/30, sigma = 1/14) # R = 1.0
param_0.9 <- c(gamma = 1/30, beta = 0.9/30, sigma = 1/14) # R = 0.9
param_0.8 <- c(gamma = 1/30, beta = 0.8/30, sigma = 1/14) # R = 0.8
# N1 <- ? # 제주도 인구
# N2 <- ? # 3주 체류 입도 인원
N <- N1 + N2
# p6 <- ? # 6월 22일 18시 기준 백신 1차 접종 인원
p <- 0.8*(p6/N)
I1 <- init_i # 7월 이동평균의 마지막 값. 6월 17일!
I2 <- 0
# r <- ? # 6월 23일자 누적격리해제 및 사망자 수
e <- 0
init_state_values <- c(S1 = (N1-I1-e-r)*(1-p),
                       S2 = N2,
                       E = e,
                       I = I1+I2,
                       R = (N1-I1-e-r)*p+r)

## (3) 초기 예측값 생성
make_initial <- function(.parameters){
    ode(y = init_state_values, func = jeju_seir, parms = .parameters, time = 1:7) %>% # 일주일 간격 갱신
        as_tibble %>% 
        pivot_longer(S1:R, names_to = "State") %>% 
        arrange(desc(State)) %>% 
        mutate(value = as.numeric(value),
               time = as.numeric(time))
}

seir1.3_init <- make_initial(param_1.3)
seir1.2_init <- make_initial(param_1.2)
seir1.1_init <- make_initial(param_1.1)
seir1_init <- make_initial(param_1.0)
seir0.9_init <- make_initial(param_0.9)
seir0.8_init <- make_initial(param_0.8)
seir_all_init <- list(seir1.3_init, seir1.2_init, seir1.1_init, 
                      seir1_init, seir0.9_init, seir0.8_init)

make_weekly <- function(seir_init, p0, parameters){ 
    
    p <- 0.8*(p0/N)
    S1 <- seir_init %>% 
        filter(State == "S1") %>%  
        mutate(value = as.numeric(value)) %>% 
        pull(value) %>% 
        tail(1)
    S2 <- seir_init %>% 
        filter(State == "S2") %>%  
        mutate(value = as.numeric(value)) %>% 
        pull(value) %>% 
        tail(1)
    i <- seir_init %>% 
        filter(State == "I") %>%  
        mutate(value = as.numeric(value)) %>% 
        pull(value) %>% 
        tail(1)
    e <- seir_init %>% 
        filter(State == "E") %>%  
        mutate(value = as.numeric(value)) %>% 
        pull(value) %>% 
        tail(1)
    r <- seir_init %>% 
        filter(State == "R") %>%  
        mutate(value = as.numeric(value)) %>% 
        pull(value) %>% 
        tail(1)
    init_state_values <- c(S1 = S1*(1-p),
                           S2 = S2,
                           E = e,
                           I = i,
                           R = S1*p+r)
    
    ode(y = init_state_values, func = jeju_seir, parms = parameters, times = 1:7) %>%
        as_tibble %>% 
        pivot_longer(S1:R, names_to = "State") %>% 
        arrange(desc(State)) %>% 
        mutate(value = as.numeric(value),
               time = as.numeric(time))
}

## (4) 일주일 간격 갱신
# 14000명씩. 15주. 21만명 추가접종. 그럼 10월 1일까지 40만명 접종완료.(9월 말까지 제주인구 70%인 40만명 접종완료계획)
# 약 15주 예측
per_week <- 1:15
.p0 <- p6+per_week*14000
seir <- vector(mode = "list", length = length(per_week))
.param <- list(param_1.3, param_1.2, param_1.1, param_1.0, param_0.9, param_0.8)
names(.param) <- str_c("seir", seq(1.3, 0.8, -0.1))
seir_all <- vector(mode = "list", length = length(.param))
names(seir_all) <- names(.param)

for(i in seq_along(.param)){
    
    for(j in seq_along(seir)){
        # SEIR 모형을 이용한 예측값 저장 공간
        if(j == 1){ # 초깃값
            seir[[j]] <- make_weekly(seir_init = seir_all_init[[i]], 
                                     p0 = .p0[[j]],
                                     parameters = .param[[i]])
        }else {
            seir[[j]] <- make_weekly(seir_init = seir[[j-1]], 
                                     p0 = .p0[[j]], 
                                     parameters = .param[[i]])
        }
        
    }
    
    seir_all[[i]] <- seir # 기초재생산지수별 SEIR 모형의 예측값 각 리스트객체에 저장.
    
}

make_final_data <- function(seir_5, seir_all){
    bind_rows(seir_5, list.rbind(seir_all)) %>% 
        arrange(desc(State)) %>% 
        mutate(stateDt = rep(seq(ymd("2021-06-18"), by = "days", length.out = (1+15)*7), 5))
}

seir_final <- 
    bind_rows(
        make_final_data(seir1.3_init, seir_all[[1]]),
        make_final_data(seir1.2_init, seir_all[[2]]),
        make_final_data(seir1.1_init, seir_all[[3]]),
        make_final_data(seir1_init, seir_all[[4]]),
        make_final_data(seir0.9_init, seir_all[[5]]),
        make_final_data(seir0.8_init, seir_all[[6]])
    ) %>% 
    mutate(R = factor(rep(seq(1.3, 0.8, by = -0.1), each = (1+15)*7*5))) %>% 
    filter(State == "I", stateDt < ymd("2021-10-01")) # 10월 전(9월 30일)까지 예측 수행


## (5) 시나리오 2: 7월 1일부로 방역 완화. R값이 1.1배씩 증가.
per_week <- 1:15
.p0 <- p6+per_week*14000
seir <- vector(mode = "list", length = length(per_week))
.param <- list(param_1.3, param_1.2, param_1.1, param_1.0, param_0.9, param_0.8)
names(.param) <- str_c("seir", seq(1.3, 0.8, -0.1))
seir_all <- vector(mode = "list", length = length(.param))
names(seir_all) <- names(.param)
for(i in seq_along(.param)){
    
    for(j in seq_along(seir)){
        # SEIR 모형을 이용한 예측값 저장 공간
        if(j == 1){ # 초깃값
            seir[[j]] <- make_weekly(seir_init = seir_all_init[[i]], 
                                     p0 = .p0[[j]],
                                     parameters = .param[[i]])
        }else {
            .param[[i]][[2]] <- .param[[i]][[2]]*1.1 # 7월 1일부터 방역완화: 재생산수 매주 1.1배
            seir[[j]] <- make_weekly(seir_init = seir[[j-1]], 
                                     p0 = .p0[[j]], 
                                     parameters = .param[[i]])
        }
        
    }
    
    seir_all[[i]] <- seir # 기초재생산지수별 SEIR 모형의 예측값 각 리스트객체에 저장.
    
}
seir_final_adjust <- 
    bind_rows(
        make_final_data(seir1.3_init, seir_all[[1]]),
        make_final_data(seir1.2_init, seir_all[[2]]),
        make_final_data(seir1.1_init, seir_all[[3]]),
        make_final_data(seir1_init, seir_all[[4]]),
        make_final_data(seir0.9_init, seir_all[[5]]),
        make_final_data(seir0.8_init, seir_all[[6]])
    ) %>% 
    mutate(R = factor(rep(seq(1.3, 0.8, by = -0.1), each = (1+15)*7*5))) %>% 
    filter(State == "I", stateDt < ymd("2021-10-01")) # 10월 전(9월 30일)까지 예측 수행

## (5) 시나리오 3: 7월 1일부로 방역 완화. 7월부터 갑자기 R값 2배 증가.
per_week <- 1:15
.p0 <- p6+per_week*14000
seir <- vector(mode = "list", length = length(per_week))
.param <- list(param_1.3, param_1.2, param_1.1, param_1.0, param_0.9, param_0.8)
names(.param) <- str_c("seir", seq(1.3, 0.8, -0.1))
seir_all <- vector(mode = "list", length = length(.param))
names(seir_all) <- names(.param)
for(i in seq_along(.param)){
    
    for(j in seq_along(seir)){
        # SEIR 모형을 이용한 예측값 저장 공간
        if(j == 1){ # 초깃값
            seir[[j]] <- make_weekly(seir_init = seir_all_init[[i]], 
                                     p0 = .p0[[j]],
                                     parameters = .param[[i]])
        }else if(j == 2) {
            .param[[i]][[2]] <- .param[[i]][[2]]*2
            seir[[j]] <- make_weekly(seir_init = seir[[j-1]], 
                                     p0 = .p0[[j]], 
                                     parameters = .param[[i]])
        }else{
            seir[[j]] <- make_weekly(seir_init = seir[[j-1]], 
                                     p0 = .p0[[j]], 
                                     parameters = .param[[i]])
        }
        
    }
    
    seir_all[[i]] <- seir # 기초재생산지수별 SEIR 모형의 예측값 각 리스트객체에 저장.
    
}
seir_final_adjust2 <- 
    bind_rows(
        make_final_data(seir1.3_init, seir_all[[1]]),
        make_final_data(seir1.2_init, seir_all[[2]]),
        make_final_data(seir1.1_init, seir_all[[3]]),
        make_final_data(seir1_init, seir_all[[4]]),
        make_final_data(seir0.9_init, seir_all[[5]]),
        make_final_data(seir0.8_init, seir_all[[6]])
    ) %>% 
    mutate(R = factor(rep(seq(1.3, 0.8, by = -0.1), each = (1+15)*7*5))) %>% 
    filter(State == "I", stateDt < ymd("2021-10-01")) # 10월 전(9월 30일)까지 예측 수행

## (5) 시나리오 4: 7월 1일부로 방역 완화. 7월부터 갑자기 R값 2배 증가 및 매주 1.1배
per_week <- 1:15
.p0 <- p6+per_week*14000
seir <- vector(mode = "list", length = length(per_week))
.param <- list(param_1.3, param_1.2, param_1.1, param_1.0, param_0.9, param_0.8)
names(.param) <- str_c("seir", seq(1.3, 0.8, -0.1))
seir_all <- vector(mode = "list", length = length(.param))
names(seir_all) <- names(.param)
for(i in seq_along(.param)){
    
    for(j in seq_along(seir)){
        # SEIR 모형을 이용한 예측값 저장 공간
        if(j == 1){ # 초깃값
            seir[[j]] <- make_weekly(seir_init = seir_all_init[[i]], 
                                     p0 = .p0[[j]],
                                     parameters = .param[[i]])
        }else if(j == 2) {
            .param[[i]][[2]] <- .param[[i]][[2]]*2
            seir[[j]] <- make_weekly(seir_init = seir[[j-1]], 
                                     p0 = .p0[[j]], 
                                     parameters = .param[[i]])
        }else{
            .param[[i]][[2]] <- .param[[i]][[2]]*1.1
            seir[[j]] <- make_weekly(seir_init = seir[[j-1]], 
                                     p0 = .p0[[j]], 
                                     parameters = .param[[i]])
        }
        
    }
    
    seir_all[[i]] <- seir # 기초재생산지수별 SEIR 모형의 예측값 각 리스트객체에 저장.
    
}
seir_final_adjust3 <- 
    bind_rows(
        make_final_data(seir1.3_init, seir_all[[1]]),
        make_final_data(seir1.2_init, seir_all[[2]]),
        make_final_data(seir1.1_init, seir_all[[3]]),
        make_final_data(seir1_init, seir_all[[4]]),
        make_final_data(seir0.9_init, seir_all[[5]]),
        make_final_data(seir0.8_init, seir_all[[6]])
    ) %>% 
    mutate(R = factor(rep(seq(1.3, 0.8, by = -0.1), each = (1+15)*7*5))) %>% 
    filter(State == "I", stateDt < ymd("2021-10-01")) # 10월 전(9월 30일)까지 예측 수행

# 3. Visualizing ----------------------------------------------------------
jeju2 %>%
    filter(!is.na(case_MA)) %>% 
    ggplot() +
    geom_line(aes(x = date, y = case_MA), col = "grey") +
    geom_line(data = seir_final,
              aes(x = stateDt, y = value, col = R)) +
    scale_x_date(date_breaks = "7 days",
                 date_labels = "%Y-%m-%d") +
    scale_y_continuous(breaks = seq(0, 16, 4),
                       limits = c(0, 16)) +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "bottom"
    ) +
    guides(colorbar = FALSE) +
    labs(
        x = "",
        y = "Daily confirmed cases",
        col = "Basic reproduction number",
        title = "SEIR model for Jeju Island: Scenario 1",
        subtitle = "2021-06-17~2021-09-30" # 7-MA 기준으로 plotting
    )
ggsave("./plot/SEIR_jeju.png", dpi = 130, width = 7.5,device = "png", type = "cairo")

jeju2 %>%
    filter(!is.na(case_MA)) %>% 
    ggplot() +
    geom_line(aes(x = date, y = case_MA), col = "grey") +
    geom_line(data = seir_final_adjust,
              aes(x = stateDt, y = value, col = R)) +
    geom_vline(xintercept = ymd("2021-07-01"), col = "blue", linetype = "dashed") +
    scale_x_date(date_breaks = "7 days",
                 date_labels = "%Y-%m-%d") +
    scale_y_continuous(breaks = seq(0, 16, 4),
                       limits = c(0, 16)) +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "bottom"
    ) +
    guides(colorbar = FALSE) +
    labs(
        x = "",
        y = "Daily confirmed cases",
        col = "Initial basic reproduction number",
        title = "SEIR model for Jeju Island: Scenario 2",
        subtitle = "2021-06-17~2021-09-30" # 7-MA 기준으로 plotting
    )
ggsave("./plot/SEIR_jeju_방역완화.png", dpi = 130, width = 7.5, device = "png", type = "cairo")

jeju2 %>%
    filter(!is.na(case_MA)) %>% 
    ggplot() +
    geom_line(aes(x = date, y = case_MA), col = "grey") +
    geom_line(data = seir_final_adjust2,
              aes(x = stateDt, y = value, col = R)) +
    geom_vline(xintercept = ymd("2021-07-01"), col = "blue", linetype = "dashed") +
    scale_x_date(date_breaks = "7 days",
                 date_labels = "%Y-%m-%d") +
    scale_y_continuous(breaks = seq(0, 16, 4),
                       limits = c(0, 16)) +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "bottom"
    ) +
    guides(colorbar = FALSE) +
    labs(
        x = "",
        y = "Daily confirmed cases",
        col = "Initial basic reproduction number",
        title = "SEIR model for Jeju Island: Scenario 3",
        subtitle = "2021-06-17~2021-09-30" # 7-MA 기준으로 plotting
    )
ggsave("./plot/SEIR_jeju_방역완화2.png", dpi = 130, width = 7.5, device = "png", type = "cairo")

jeju2 %>%
    filter(!is.na(case_MA)) %>% 
    ggplot() +
    geom_line(aes(x = date, y = case_MA), col = "grey") +
    geom_line(data = seir_final_adjust3,
              aes(x = stateDt, y = value, col = R)) +
    geom_vline(xintercept = ymd("2021-07-01"), col = "blue", linetype = "dashed") +
    scale_x_date(date_breaks = "7 days",
                 date_labels = "%Y-%m-%d") +
    scale_y_continuous(breaks = seq(0, 16, 4),
                       limits = c(0, 16)) +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "bottom"
    ) +
    guides(colorbar = FALSE) +
    labs(
        x = "",
        y = "Daily confirmed cases",
        col = "Initial basic reproduction number",
        title = "SEIR model for Jeju Island: Scenario 4",
        subtitle = "Weekly update basic reproduction number from Jul 1, 2021\n2021-06-17~2021-09-30" # 7-MA 기준으로 plotting
    )
ggsave("./plot/SEIR_jeju_방역완화3_specification.png", dpi = 130, width = 7.5, device = "png", type = "cairo")
