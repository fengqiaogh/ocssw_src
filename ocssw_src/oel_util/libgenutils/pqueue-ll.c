#include "pqueue.h"

#include <stdio.h>
#include <stdlib.h>

typedef struct pqueue_node pqueue_node;

struct pqueue_node {
    double priority;
    void *item;
    pqueue_node *next;
};

struct pqueue {
    pqueue_node *start;
};

pqueue *pqueue_create() {
    pqueue *ret = (pqueue*) malloc(sizeof(pqueue));
    ret->start = NULL;
    return ret;
}

void pqueue_destroy(pqueue *pq) {
    pqueue_node *tmp = pq->start;
    while (tmp) {
        pqueue_node *tmp_was = tmp;
        tmp = tmp->next;
        free(tmp_was);
    };
    free(pq);
}

void pqueue_push(pqueue *pq, double priority, void *item) {
    pqueue_node *new = (pqueue_node*) malloc(sizeof(pqueue_node));
    new->priority = priority;
    new->item = item;
    new->next = NULL;
    if (pq->start == NULL) {
        pq->start = new;
    } else if (priority <= pq->start->priority) {
        new->next = pq->start;
        pq->start = new;
    } else {
        pqueue_node *q = pq->start;
        while (q->next != NULL && q->next->priority <= priority) {
            q = q->next;
        }
        new->next = q->next;
        q->next = new;
    }
}

void pqueue_repush(pqueue *pq, double priority, void *item) {
    pqueue_node *current, *current_parent = NULL;
    for (current = pq->start; current != NULL; current_parent = current, current = current->next) {
        if (current->item == item) {
            if (current->priority != priority) {
                if (pq->start->priority > priority) {
                    current->priority = priority;
                    if (current_parent) {
                        current_parent->next = current->next;
                        current->next = pq->start;
                        pq->start = current;
                    }
                } else {
                    pqueue_node *new_parent;
                    if (priority > current->priority) {
                        new_parent = current;
                    } else {
                        new_parent = pq->start;
                    }
                    while (new_parent->next != NULL && new_parent->next->priority <= priority) {
                        new_parent = new_parent->next;
                    }
                    current->priority = priority;
                    if (new_parent != current) {
                        if (current_parent) {
                            current_parent->next = current->next;
                        } else {
                            pq->start = pq->start->next;
                        }
                        pqueue_node *tmp = new_parent->next;
                        new_parent->next = current;
                        current->next = tmp;
                    }
                }
                break;
            }
        }
    }
}

void* pqueue_pull(pqueue *pq) {
    if (pq->start == NULL) {
        return NULL;
    } else {
        pqueue_node *old_start = pq->start;
        pq->start = pq->start->next;
        void *tmp = old_start->item;
        free(old_start);
        return tmp;
    }
}

void pqueue_display(pqueue *pq, void (*print_item)(void *item)) {
    printf("[");
    pqueue_node *tmp;
    for (tmp = pq->start; tmp != NULL; tmp = tmp->next) {
        printf("p%.2f => ", tmp->priority);
        if (print_item) {
            print_item(tmp->item);
        } else {
            printf("%p", tmp->item);
        }
        if (tmp->next != NULL) {
            printf(", ");
        }
    }
    printf("]");
}
